# c6_coherence.R -- fewclust 0.3.0, prediction P1 (pre-registered in dev/aggiornamento_0.3.0.md): how often the
# coherence condition of fewclust_control(check_coherence = TRUE), which is not part of the paper's rule, changes the
# reported marginalized-model fit on study datasets. Pool B: scenario M3 (latent standard deviation growing over time),
# four mechanisms, N in {25, 65}, T in {4, 8}, H0, replicates 1-50 (800 datasets), regenerated from their seeds with
# the code of branch morris-v4 (read only). Writes only to <out_dir>.
# Usage: Rscript --vanilla dev/checks_0.3.0/c6_coherence.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROJ <- "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
.libPaths(c(DEV, PROJ)); invisible(loadNamespace("fewclust")); stopifnot(getNamespaceVersion("fewclust") == "0.3.0")
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))  # generators and seeds
# the study code calls library(fewclust), which under R 4.6 replaces the loaded 0.3.0 with the 0.2.0 of ~/R/s7d_lib:
# detach it and load the development version again (the study functions used here do not call fewclust)
if ("package:fewclust" %in% search()) detach("package:fewclust", unload = TRUE)
if ("fewclust" %in% loadedNamespaces()) unloadNamespace("fewclust")
invisible(loadNamespace("fewclust", lib.loc = DEV))
stopifnot(getNamespaceVersion("fewclust") == "0.3.0", dirname(getNamespaceInfo("fewclust", "path")) == DEV)
suppressMessages(library(data.table)); setDTthreads(1L)
TA <- s7d2_tasks()
MECH <- c("GenBinaryY", "GLMM-native", "GLMM-AR1", "Copula-AR1")
rows <- TA[dgp %in% MECH & scen == "M3" & N %in% c(25L, 65L) & ntime %in% c(4L, 8L) & label == "manual" & chunk == 1L]
stopifnot(nrow(rows) == 16L)
POOL <- rows[, .(sim_id = 1:50), by = .(dgp, scen, N, ntime, label, model, old_file, file, chunk)]
one <- function(k) {
  t <- POOL[k]; i <- t$sim_id
  gen <- s7d2_gen(t)
  invisible(design_degenerate(t$model, t$N, t$ntime, i))
  set_rep_seed(t$model, i)
  d <- as.data.frame(gen(t$N, t$ntime, i, t$label)); d <- d[order(d$id, d$time), ]
  p <- fewclust:::fc_parse(binY ~ time * trt, d, quote(id), quote(time))
  df <- nlevels(p$id) - ncol(p$X)
  fit <- function(coh) fewclust:::fit_mmm(p, fewclust::fewclust_control(check_coherence = coh), "t", 0.05, df)
  a <- fit(FALSE); b <- fit(TRUE)
  ok <- !inherits(a, "fc_failure") && !inherits(b, "fc_failure")
  ja <- if (ok) match("time:trt", a$coefficients$term) else NA
  data.table(dgp = t$dgp, N = t$N, ntime = t$ntime, sim_id = i, both_fitted = ok,
             path_published = if (ok) a$model$check$path else NA_character_,
             path_coherence = if (ok) b$model$check$path else NA_character_,
             changed = if (ok) !identical(a$coefficients, b$coefficients) else NA,
             est_published = if (ok) a$coefficients$estimate[ja] else NA_real_,
             est_coherence = if (ok) b$coefficients$estimate[ja] else NA_real_,
             se_published = if (ok) a$coefficients$std.error[ja] else NA_real_,
             se_coherence = if (ok) b$coefficients$std.error[ja] else NA_real_)
}
t0 <- proc.time()[[3]]
R <- rbindlist(parallel::mclapply(seq_len(nrow(POOL)), one, mc.cores = 4L, mc.preschedule = FALSE))
fwrite(R, file.path(OUT, "c6_coherence.csv"))
sink(file.path(OUT, "c6_summary.txt"), split = TRUE)
cat("seconds:", round(proc.time()[[3]] - t0), "; datasets:", nrow(R), "; both fits returned:", R[, sum(both_fitted)], "\n")
n <- R[both_fitted == TRUE, .N]; k <- R[both_fitted == TRUE, sum(changed)]; p <- k / n; m <- sqrt(p * (1 - p) / n)
cat(sprintf("P1: share of fits changed by the coherence condition = %d / %d = %.2f%% (MCSE %.2f; estimate +/- 2 MCSE: %.2f%% to %.2f%%)\n",
            k, n, 100 * p, 100 * m, 100 * (p - 2 * m), 100 * (p + 2 * m)))
cat("relevance interval 1% to 6%:", if (p - 2 * m >= 0.01 && p + 2 * m <= 0.06) "CONFIRMED" else
    if (p + 2 * m < 0.01 || p - 2 * m > 0.06) "REFUTED" else "UNDETERMINED", "\n")
print(R[both_fitted == TRUE, .(n = .N, changed = sum(changed)), by = .(dgp, N, ntime)])
print(R[both_fitted == TRUE, .N, by = .(path_published, path_coherence)])
ch <- R[changed == TRUE]
if (nrow(ch)) { cat("among changed fits, time:trt estimate difference / SE (published):\n")
  print(summary(ch[, (est_coherence - est_published) / se_published]))
  cat("SE ratio coherence / published:\n"); print(summary(ch[, se_coherence / se_published])) }
sink()
