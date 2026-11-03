# c9b_diagnose.R -- fewclust 0.3.1, explanation of check 3 (NOT a gate; written after check 3 failed, 2026-11-02).
# Reads out/c9_mmm_repro.csv and (1) applies the usability filter of the paper (Section 2.5: finite estimate and SE,
# |estimate| < 10, 0 < SE < 10) to the study's and the package's fit of time:trt; (2) recounts the description of the
# default control with missing values handled (the c9 summary printed NA); (3) for each mismatch, regenerates the
# dataset and records the study's fit A/B and the package's trend and constant fits (code, slope, log-likelihood, SE).
# Usage: Rscript --vanilla dev/checks_0.3.1/c9b_diagnose.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
PROJ <- "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
.libPaths(c(DEV, PROJ))
invisible(loadNamespace("fewclust"))
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))
if ("package:fewclust" %in% search()) detach("package:fewclust", unload = TRUE)
if ("fewclust" %in% loadedNamespaces()) unloadNamespace("fewclust")
invisible(loadNamespace("fewclust", lib.loc = DEV)); stopifnot(getNamespaceVersion("fewclust") == "0.3.1")
suppressMessages(library(data.table)); setDTthreads(1L)
R <- fread(file.path(OUT, "c9_mmm_repro.csv"))
usable <- function(e, s) is.finite(e) & is.finite(s) & abs(e) < 10 & s > 0 & s < 10
R[, study_usable := compared & usable(est_saved, se_saved)]
R[, pkg_usable := compared & usable(est_pkg, se_pkg)]
sink(file.path(OUT, "c9b_summary.txt"), split = TRUE)
cat("compared:", R[, sum(compared)], "; mismatches:", R[compared == TRUE, sum(!match)], "\n")
cat("usable by the paper's filter: study", R[, sum(study_usable)], "; package", R[, sum(pkg_usable)],
    "; both", R[, sum(study_usable & pkg_usable)], "; neither (among compared)", R[compared == TRUE, sum(!study_usable & !pkg_usable)], "\n")
cat("mismatches among datasets usable for the study:", R[study_usable == TRUE, sum(!match)], "of", R[, sum(study_usable)],
    "; among datasets usable for neither:", R[compared == TRUE & !study_usable & !pkg_usable, sum(!match)], "\n")
cat("usable for one and not the other:", R[compared == TRUE, sum(study_usable != pkg_usable)], "\n")
dd <- R[pkg_status == "fitted" & default_status == "fitted"]
dd[, differs := !(mapply(identical, est_pkg, est_default) & mapply(identical, se_pkg, se_default))]
cat("description, default control (check_coherence = TRUE): reported fit different in", dd[, sum(differs)], "of", nrow(dd),
    "fitted datasets; path changed in", dd[, sum(path_pkg != default_path, na.rm = TRUE)], "\n")
if (dd[, sum(differs)]) print(dd[differs == TRUE, .(dgp, N, ntime, label, sim_id, path_pkg, default_path, est_pkg, est_default, se_pkg, se_default)])
M <- R[compared == TRUE & match == FALSE]
TA <- s7d2_tasks()
det <- rbindlist(lapply(seq_len(nrow(M)), function(k) {
  m <- M[k]; t <- TA[dgp == m$dgp & scen == "M3" & N == m$N & ntime == m$ntime & label == m$label & chunk == 1L]
  cc <- readRDS(file.path("/home/user/MMM/results_v4/S6cc/rows", sprintf("%s_m3_N%d_T%d_%s_c01.rds", gsub("[^A-Za-z0-9]", "", m$dgp), m$N, m$ntime, m$label)))[sim_id == m$sim_id]
  gen <- s7d2_gen(t); invisible(design_degenerate(t$model, t$N, t$ntime, m$sim_id)); set_rep_seed(t$model, m$sim_id)
  d <- as.data.frame(gen(t$N, t$ntime, m$sim_id, t$label)); d <- d[order(d$id, d$time), ]
  f <- fewclust::auto_cluster_fit(binY ~ time * trt, d, group = "id", time = "time", control = fewclust::fewclust_control(check_coherence = FALSE))
  ck <- f$model$check
  data.table(dgp = m$dgp, N = m$N, ntime = m$ntime, label = m$label, sim_id = m$sim_id, events = sum(d$binY), n_obs = nrow(d),
             study_used = cc$used, study_A_code = cc$A_code, study_A_slope = cc$A_slope_min, study_A_est = cc$A_est, study_A_se = cc$A_se,
             study_B_code = cc$B_code, study_B_est = cc$B_est, study_B_se = cc$B_se,
             pkg_path = ck$path, pkg_trend_code = ck$trend$code, pkg_trend_slope = ck$trend$slope_min,
             pkg_trend_est = unname(ck$trend$estimate["time:trt"]), pkg_trend_se = unname(ck$trend$std.error["time:trt"]),
             pkg_const_code = if (is.null(ck$constant)) NA_integer_ else ck$constant$code,
             pkg_const_est = if (is.null(ck$constant)) NA_real_ else unname(ck$constant$estimate["time:trt"]),
             pkg_const_se = if (is.null(ck$constant)) NA_real_ else unname(ck$constant$std.error["time:trt"]),
             pkg_trend_logLik = ck$trend$logLik, pkg_const_logLik = if (is.null(ck$constant)) NA_real_ else ck$constant$logLik,
             study_usable = M$study_usable[k], pkg_usable = M$pkg_usable[k])
}))
fwrite(det, file.path(OUT, "c9b_mismatches.csv"))
options(width = 200); print(det)
sink()
