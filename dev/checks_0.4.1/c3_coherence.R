# c3_coherence.R -- fewclust 0.4.1, prediction 4 (dev/aggiornamento_0.4.1.md): the comparison of check 3 of 0.3.1
# (dev/checks_0.3.1/c9_mmm_repro.R + c9b_diagnose.R, whose counts "127 compared" and "reported fit different in 4"
# are cited by ~/MMM/R/s9/s9_fewclust_coherence.R -> tables/V6_S9_fewclust_coherence.csv) rerun with fewclust 0.4.1:
# same pool (scenario M3, GenBinaryY, GLMM-native, GLMM-AR1, N in {10, 25, 65}, T in {2, 4, 8}, H0 and H1, replicates
# 1-3: 162 datasets) regenerated with the study's code (read only), same study fits (results_v4/S6cc/rows), same
# definitions; the call auto_cluster_fit(binY ~ time * trt, d, group = "id", time = "time", control = ctrl) of 0.3.1
# becomes fit_fewclust(d, binY ~ time * trt, id = "id", time = "time", control = ctrl) (method "auto").
# Items with suffix _d90: the same counts against the study fit under D-90 (study gives no result when fit A is at the
# boundary and fit B is not finite or has nlm code not 1-2).
# Usage (from the package root): Rscript --vanilla dev/checks_0.4.1/c3_coherence.R <lib> <out_dir>
args <- commandArgs(TRUE); FC_LIB <- normalizePath(args[1]); OUT <- args[2]
PROJ <- "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
.libPaths(c(FC_LIB, PROJ))
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))
# the study code also resets .libPaths() and defines globals (LIB among them): use FC_LIB and set the paths again
.libPaths(c(FC_LIB, PROJ))
if ("package:fewclust" %in% search()) detach("package:fewclust", unload = TRUE)
if ("fewclust" %in% loadedNamespaces()) unloadNamespace("fewclust")
invisible(loadNamespace("fewclust", lib.loc = FC_LIB))
stopifnot(getNamespaceVersion("fewclust") == "0.4.1", dirname(getNamespaceInfo("fewclust", "path")) == FC_LIB)
suppressMessages(library(data.table)); setDTthreads(1L)
fc <- asNamespace("fewclust")
TA <- s7d2_tasks()
MECH <- c("GenBinaryY", "GLMM-native", "GLMM-AR1")
rows <- TA[dgp %in% MECH & scen == "M3" & N %in% c(10L, 25L, 65L) & ntime %in% c(2L, 4L, 8L) & chunk == 1L]
stopifnot(nrow(rows) == 54L)
CC <- "/home/user/MMM/results_v4/S6cc/rows"
rows[, cc_file := file.path(CC, sprintf("%s_m%d_N%d_T%d_%s_c01.rds", gsub("[^A-Za-z0-9]", "", dgp), model, N, ntime, label))]
stopifnot(all(file.exists(rows$cc_file)))
POOL <- rows[rep(seq_len(.N), each = 3L)][, sim_id := rep(1:3, .N / 3L)]
files <- unique(c(POOL$file, POOL$cc_file))
md5_before <- tools::md5sum(files)
paper <- fc$fewclust_control(check_coherence = FALSE)
one <- function(k) {
  t <- POOL[k]; i <- t$sim_id
  gen <- s7d2_gen(t)
  saved <- readRDS(t$file)[sim_id == i]
  cc <- readRDS(t$cc_file)[sim_id == i]
  stopifnot(nrow(saved) == 1L, nrow(cc) <= 1L)
  base <- data.table(dgp = t$dgp, N = t$N, ntime = t$ntime, label = t$label, sim_id = i)
  deg <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) e)
  if (inherits(d, "error")) return(cbind(base, design_degenerate = deg, generated = FALSE))
  d <- d[order(d$id, d$time), ]
  cc_ok <- nrow(cc) == 1L && isTRUE(cc$A_ok) && !isFALSE(cc$A_has_term) && is.finite(cc$cc_est) && is.finite(cc$cc_se)
  study_bnd <- nrow(cc) == 1L && isTRUE(cc$bnd)
  study_Bconv <- nrow(cc) == 1L && isTRUE(cc$B_ok) && is.finite(cc$B_est) && is.finite(cc$B_se) && cc$B_code %in% 1:2
  fit_of <- function(ctrl) tryCatch(suppressWarnings(suppressMessages(
    fc$fit_fewclust(d, binY ~ time * trt, id = "id", time = "time", control = ctrl))), error = function(e) e)
  a <- fit_of(paper); b <- fit_of(fc$fewclust_control())
  get <- function(f, what) {
    if (!inherits(f, "fewclust")) return(NA)
    tab <- f$table; j <- match("time:trt", tab$term)
    switch(what, strategy = f$method, est = tab$estimate[j], se = tab$std.error[j],
           path = if (!is.null(f$convergence)) f$convergence$path else NA_character_)
  }
  status <- function(f) if (inherits(f, "fewclust")) "fitted" else paste("error:", conditionMessage(f))
  cbind(base, design_degenerate = deg, generated = TRUE,
        data_identical = nrow(d) == saved$n_obs && sum(d$binY) == saved$sum_y,
        constant = sum(d$binY) %in% c(0, nrow(d)),
        cc_ok = cc_ok, study_no_result_d90 = study_bnd && !study_Bconv,
        used_saved = if (nrow(cc)) cc$used else NA_character_,
        code_saved = if (!nrow(cc) || is.na(cc$used)) NA_integer_ else if (cc$used == "B") cc$B_code else cc$A_code,
        est_saved = if (nrow(cc)) cc$cc_est else NA_real_, se_saved = if (nrow(cc)) cc$cc_se else NA_real_,
        pkg_status = status(a), pkg_strategy = get(a, "strategy"), path_pkg = get(a, "path"),
        est_pkg = get(a, "est"), se_pkg = get(a, "se"),
        default_status = status(b), default_strategy = get(b, "strategy"), default_path = get(b, "path"),
        est_default = get(b, "est"), se_default = get(b, "se"))
}
one_safe <- function(k) tryCatch(one(k), error = function(e)
  data.table(dgp = POOL$dgp[k], N = POOL$N[k], ntime = POOL$ntime[k], label = POOL$label[k], sim_id = POOL$sim_id[k],
             error = conditionMessage(e)))
t0 <- proc.time()[[3]]
R <- rbindlist(parallel::mclapply(seq_len(nrow(POOL)), one_safe, mc.cores = 5L, mc.preschedule = FALSE), fill = TRUE)
if (!"error" %in% names(R)) R[, error := NA_character_]
secs <- round(proc.time()[[3]] - t0)
stopifnot(nrow(R) == nrow(POOL))
eq <- function(a, b, tol) mapply(function(x, y) isTRUE(all.equal(x, y, tolerance = tol)), a, b)
usable <- function(e, s) is.finite(e) & is.finite(s) & abs(e) < 10 & s > 0 & s < 10
# definitions of c9_mmm_repro.R / c9b_diagnose.R
R[, compared := generated %in% TRUE & cc_ok %in% TRUE & pkg_status %in% "fitted" & is.na(error)]
R[, numeric_compared := compared & code_saved %in% c(1L, 2L)]
R[, branch_match := (used_saved == "B") == (path_pkg %in% "constant fit used")]
R[, match := ifelse(!compared, NA, pkg_strategy %in% "mmm" & branch_match &
                      (!numeric_compared | (eq(est_saved, est_pkg, 1e-6) & eq(se_saved, se_pkg, 1e-5))))]
R[, study_usable := compared & usable(est_saved, se_saved)]
R[, pkg_usable := compared & usable(est_pkg, se_pkg)]
R[, both_fitted := pkg_status %in% "fitted" & default_status %in% "fitted"]
R[, differs := both_fitted & !(mapply(identical, est_pkg, est_default) & mapply(identical, se_pkg, se_default))]
# the same against the study under D-90: where the study gives no result the package must give no result of the
# marginalized model (GEE reported); elsewhere as above
R[, compared_d90 := compared & !study_no_result_d90]
R[, match_d90 := ifelse(!compared_d90, NA, match)]
R[, key_no_mmm := compared & study_no_result_d90 & !pkg_strategy %in% "mmm"]
md5_after <- tools::md5sum(files)
fwrite(R, file.path(OUT, "c3_coherence_rows.csv"))
src <- "fewclust 0.4.1: fewclust/dev/checks_0.4.1/c3_coherence.R -> out/c3_coherence_rows.csv"
X <- data.table(
  item = c("n_cmp", "n_mismatch", "n_usable_study", "n_mismatch_usable", "n_coh",
           "n_cmp_d90", "n_mismatch_d90", "n_study_no_result_d90", "n_study_no_result_d90_pkg_no_mmm"),
  value = c(R[, sum(compared)], R[compared == TRUE, sum(!match)], R[, sum(study_usable)],
            R[study_usable == TRUE, sum(!match)], R[, sum(differs)],
            R[, sum(compared_d90)], R[compared_d90 == TRUE, sum(!match_d90)], R[compared == TRUE, sum(study_no_result_d90)],
            R[, sum(key_no_mmm)]),
  source = src)
fwrite(X, file.path(OUT, "c3_coherence.csv"))
sink(file.path(OUT, "c3_summary.txt"), split = TRUE)
cat("fewclust", as.character(getNamespaceVersion("fewclust")), "; seconds:", secs, "\n")
cat("datasets:", nrow(R), "; errors in the check:", R[!is.na(error), .N], "; generated:", R[generated %in% TRUE, .N],
    "; data identical to the saved run (n_obs, sum_y):", R[data_identical %in% TRUE, .N], "\n")
cat("compared:", R[, sum(compared)], "; mismatches:", R[compared == TRUE, sum(!match)], "\n")
cat("mismatches among datasets usable for the study:", R[study_usable == TRUE, sum(!match)], "of", R[, sum(study_usable)], "\n")
cat("description, default control (check_coherence = TRUE): reported fit different in", R[, sum(differs)], "of",
    R[, sum(both_fitted)], "fitted datasets; path changed in", R[both_fitted == TRUE, sum(path_pkg != default_path, na.rm = TRUE)], "\n")
options(width = 220)
if (R[, sum(differs)]) print(R[differs == TRUE, .(dgp, N, ntime, label, sim_id, pkg_strategy, path_pkg, default_strategy,
                                                  default_path, est_pkg, est_default, se_pkg, se_default)])
cat("D-90: study without result among compared:", R[compared == TRUE, sum(study_no_result_d90)],
    "; of which the package reports no marginalized model (GEE):", R[, sum(key_no_mmm)], "\n")
cat("D-90: compared", R[, sum(compared_d90)], "; mismatches", R[compared_d90 == TRUE, sum(!match_d90)], "\n")
if (R[compared == TRUE & match == FALSE, .N]) { cat("mismatches:\n")
  print(R[compared == TRUE & match == FALSE, .(dgp, N, ntime, label, sim_id, used_saved, code_saved, study_no_result_d90,
                                                pkg_strategy, path_pkg, est_saved, est_pkg, se_saved, se_pkg)]) }
cat("study rows read unchanged (md5 before = after):", identical(md5_before, md5_after), "\n")
print(X)
sink()
