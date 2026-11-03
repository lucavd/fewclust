# c2_study_sample.R -- fewclust 0.4.1, check C (dev/aggiornamento_0.4.1.md): stratified sample of study datasets of
# scenarios M3-M5. Keys: results_v4/S11/cv_keys.rds (the marginalized model with the convergence check gives no result
# under D-90), strata mechanism x scenario x why, 30 per stratum in M3 and 10 in M4-M5. Converged refits: fit A at the
# boundary, fit B finite with nlm code 1-2 (used = "B"), not a key; strata mechanism x scenario, 60 in M3, 20 in M4-M5.
# Study fits: results_v4/S6cc/cc_all.rds (GenBinaryY, GLMM-native, GLMM-AR1) and results_v4/S7e/rows/*.rds $cc
# (Copula-AR1), read only, with the definitions of R/s11/s11_cv_keys.R. Seed 20261103.
# Usage: Rscript --vanilla c2_study_sample.R <out_dir> <scratch_dir>
args <- commandArgs(TRUE); OUT <- args[1]; SCR <- args[2]
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
ST <- "/home/user/MMM/results_v4"
K6 <- c("dgp", "scen", "N", "ntime", "label", "sim_id")
f_keys <- file.path(ST, "S11/cv_keys.rds"); f_cc <- file.path(ST, "S6cc/cc_all.rds")
f_cop <- list.files(file.path(ST, "S7e/rows"), "^copula_m[345]_.*\\.rds$", full.names = TRUE)
stopifnot(length(f_cop) == 3L * 16L * 6L * 2L * 10L)
md5_before <- tools::md5sum(c(f_keys, f_cc, f_cop))
KEYS <- as.data.table(readRDS(f_keys))[scen %in% c("M3", "M4", "M5")]
cols <- c(K6, "A_ok", "A_has_term", "A_code", "A_slope_min", "A_est", "A_se", "bnd", "B_ok", "B_est", "B_se", "B_code",
          "used", "cc_est", "cc_se")
A <- as.data.table(readRDS(f_cc))[scen %in% c("M3", "M4", "M5"), ..cols]
C <- rbindlist(parallel::mclapply(f_cop, function(f) as.data.table(readRDS(f)$cc), mc.cores = 5L), fill = TRUE)
C <- C[, ..cols]
X <- rbind(A, C); rm(A, C)
stopifnot(!anyDuplicated(X, by = K6))
X[, fitA := A_ok %in% TRUE & A_has_term %in% TRUE]
X[, bnd2 := fitA & (!(A_code %in% 1:2) | (is.finite(A_slope_min) & A_slope_min < 1e-3))]
X[, Bfin := B_ok %in% TRUE & is.finite(B_est) & is.finite(B_se)]
X[, Bconv := Bfin & B_code %in% 1:2]
X[, key := bnd2 & !Bconv]
# the keys file is the set of these datasets (check, not assumed)
chk <- merge(KEYS[, c(K6, "why"), with = FALSE], X[key == TRUE, c(K6, "Bfin"), with = FALSE], by = K6, all = TRUE)
stopifnot(nrow(chk) == nrow(KEYS), !anyNA(chk$why), !anyNA(chk$Bfin),
          all(chk$why == fifelse(chk$Bfin, "refit not converged", "refit not finite")))
set.seed(20261103L)
setkeyv(X, K6)
take <- function(D, by, n_m3, n_other)
  D[D[, .I[sort(sample(.N, min(.N, if (scen[1] == "M3") n_m3 else n_other)))], by = by]$V1]
KX <- merge(X[key == TRUE], KEYS[, c(K6, "why", "src"), with = FALSE], by = K6)
setkeyv(KX, K6)
SK <- take(KX, c("dgp", "scen", "why"), 30L, 10L)[, kind := "key"]
CX <- X[bnd2 == TRUE & Bconv == TRUE & used == "B" & key == FALSE]
SC <- take(CX, c("dgp", "scen"), 60L, 20L)[, `:=`(kind = "conv", why = "refit converged",
                                                     src = fifelse(dgp == "Copula-AR1", "copula", "main"))]
nm <- c(K6, "kind", "why", "src", setdiff(cols, K6))
SAMPLE <- rbind(SK[, ..nm], SC[, ..nm])
setorderv(SAMPLE, c("kind", K6))
fwrite(SAMPLE, file.path(OUT, "c2_sample.csv"))
md5_after <- tools::md5sum(c(f_keys, f_cc, f_cop))
fwrite(data.table(file = names(md5_before), md5 = unname(md5_before)), file.path(SCR, "c2_md5_sources_before.csv"))
sink(file.path(OUT, "c2_sample_summary.txt"), split = TRUE)
cat("keys in M3-M5:", nrow(KEYS), "; converged-refit boundary datasets (not keys) in M3-M5:", nrow(CX), "\n")
print(SAMPLE[, .N, by = .(kind, scen, dgp, why)][order(kind, scen, dgp, why)])
cat("sample:", nrow(SAMPLE), "datasets; keys", SAMPLE[kind == "key", .N], "; converged", SAMPLE[kind == "conv", .N],
    "; of which M3:", SAMPLE[scen == "M3" & kind == "key", .N], "and", SAMPLE[scen == "M3" & kind == "conv", .N], "\n")
cat("source files unchanged while sampling (md5):", identical(md5_before, md5_after), "\n")
sink()
