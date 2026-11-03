# c9_mmm_repro.R -- fewclust 0.3.1, check 3 (pre-registered in dev/aggiornamento_0.3.1.md): the recommended strategy of
# 0.3.1 (marginalized multilevel model with the convergence check, "mmm") reproduces the study's stored fit of the
# marginalized model with the convergence check (results_v4/S6cc/rows on branch morris-v4: MMLB::mm, q = 20, rule of
# D-32; cc_est, cc_se, used). Pool: scenario M3 (mean binY ~ time * trt, latent SD ~ time, transition ~ 1, the model the
# package fits), mechanisms GenBinaryY, GLMM-native, GLMM-AR1, N in {10, 25, 65}, T in {2, 4, 8}, H0 and H1, replicates
# 1-3 (162 datasets), regenerated from their seeds with the code of morris-v4 (read only). Data identity is checked
# against n_obs and sum_y of results_v4/S7d2/rows. Writes only to <out_dir>.
# Usage (from the package root): Rscript --vanilla dev/checks_0.3.1/c9_mmm_repro.R <dev_lib> <out_dir>
args <- commandArgs(TRUE)
DEV <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROJ <- "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
.libPaths(c(DEV, PROJ))
invisible(loadNamespace("fewclust")); stopifnot(getNamespaceVersion("fewclust") == "0.3.1")
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))  # generators, seeds, tasks
# the study code calls library(fewclust), which under R 4.6 replaces the loaded 0.3.1 with the 0.2.0 of ~/R/s7d_lib:
# detach it and load the development version again (the study functions used here do not call fewclust)
if ("package:fewclust" %in% search()) detach("package:fewclust", unload = TRUE)
if ("fewclust" %in% loadedNamespaces()) unloadNamespace("fewclust")
invisible(loadNamespace("fewclust", lib.loc = DEV))
stopifnot(getNamespaceVersion("fewclust") == "0.3.1", dirname(getNamespaceInfo("fewclust", "path")) == DEV)
suppressMessages(library(data.table)); setDTthreads(1L)

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
paper <- fewclust::fewclust_control(check_coherence = FALSE)

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
  if (inherits(d, "error")) return(cbind(base, design_degenerate = deg, generated = FALSE, generated_saved = saved$generated))
  d <- d[order(d$id, d$time), ]
  cc_ok <- nrow(cc) == 1L && isTRUE(cc$A_ok) && !isFALSE(cc$A_has_term) && is.finite(cc$cc_est) && is.finite(cc$cc_se)
  fit_of <- function(ctrl) tryCatch(fewclust::auto_cluster_fit(binY ~ time * trt, d, group = "id", time = "time", control = ctrl),
                                    error = function(e) e)
  a <- fit_of(paper); b <- fit_of(fewclust::fewclust_control())
  get <- function(f, what) {
    if (!inherits(f, "fewclust_fit")) return(NA)
    j <- match("time:trt", f$coefficients$term)
    switch(what, strategy = f$strategy, est = f$coefficients$estimate[j], se = f$coefficients$std.error[j],
           path = if (f$strategy == "mmm") f$model$check$path else NA_character_)
  }
  status <- function(f) if (inherits(f, "fewclust_fit")) "fitted" else if (inherits(f, "fewclust_refusal"))
    paste("refused:", f$decision$rule) else paste("error:", conditionMessage(f))
  cbind(base, design_degenerate = deg, generated = TRUE, generated_saved = saved$generated,
        data_identical = nrow(d) == saved$n_obs && sum(d$binY) == saved$sum_y,
        constant = sum(d$binY) %in% c(0, nrow(d)),
        cc_ok = cc_ok, used_saved = if (nrow(cc)) cc$used else NA_character_,
        code_saved = if (!nrow(cc) || is.na(cc$used)) NA_integer_ else if (cc$used == "B") cc$B_code else cc$A_code,
        est_saved = if (nrow(cc)) cc$cc_est else NA_real_, se_saved = if (nrow(cc)) cc$cc_se else NA_real_,
        pkg_status = status(a), pkg_strategy = get(a, "strategy"), path_pkg = get(a, "path"),
        est_pkg = get(a, "est"), se_pkg = get(a, "se"),
        default_status = status(b), default_path = get(b, "path"), est_default = get(b, "est"), se_default = get(b, "se"))
}
one_safe <- function(k) tryCatch(one(k), error = function(e)
  data.table(dgp = POOL$dgp[k], N = POOL$N[k], ntime = POOL$ntime[k], label = POOL$label[k], sim_id = POOL$sim_id[k],
             error = conditionMessage(e)))
t0 <- proc.time()[[3]]
R <- rbindlist(parallel::mclapply(seq_len(nrow(POOL)), one_safe, mc.cores = 4L, mc.preschedule = FALSE), fill = TRUE)
if (!"error" %in% names(R)) R[, error := NA_character_]
secs <- round(proc.time()[[3]] - t0)
stopifnot(nrow(R) == nrow(POOL))
R[, compared := generated %in% TRUE & cc_ok %in% TRUE & pkg_status %in% "fitted" & is.na(error)]
R[, numeric_compared := compared & code_saved %in% c(1L, 2L)]
R[, branch_match := (used_saved == "B") == (path_pkg %in% "constant fit used")]
R[, match := ifelse(!compared, NA, pkg_strategy %in% "mmm" & branch_match &
                      (!numeric_compared | (mapply(function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-6)), est_saved, est_pkg) &
                                            mapply(function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-5)), se_saved, se_pkg))))]
R[, default_differs := pkg_status %in% "fitted" & default_status %in% "fitted" & (est_pkg != est_default | se_pkg != se_default)]
md5_after <- tools::md5sum(files)
fwrite(R, file.path(OUT, "c9_mmm_repro.csv"))
sink(file.path(OUT, "c9_summary.txt"), split = TRUE)
cat("fewclust", as.character(getNamespaceVersion("fewclust")), "; seconds:", secs, "\n")
cat("datasets:", nrow(R), "; errors in the check:", R[!is.na(error), .N], "; generated:", R[generated %in% TRUE, .N],
    "; data identical to the saved run (n_obs, sum_y):", R[data_identical %in% TRUE, .N], "of", R[generated %in% TRUE, .N], "\n")
if (R[!is.na(error), .N]) print(R[!is.na(error), .(dgp, N, ntime, label, sim_id, error)])
cat("compared (study fit usable, package fitted):", R[, sum(compared)], "; of which numerically (study's reported fit with nlm code 1-2):",
    R[, sum(numeric_compared)], "\n")
cat("check 3 (strategy mmm, same branch, estimate all.equal 1e-6, SE all.equal 1e-5):", R[compared == TRUE, sum(match)], "of",
    R[, sum(compared)], "\n")
cat("  branch: constant fit used by the study", R[compared == TRUE, sum(used_saved == "B")], "; by the package",
    R[compared == TRUE, sum(path_pkg %in% "constant fit used")], "\n")
cat("  not compared:\n"); print(R[compared == FALSE | is.na(compared), .N, by = .(generated, cc_ok, pkg_status, constant)])
if (R[compared == TRUE & match == FALSE, .N]) { cat("  mismatches:\n")
  print(R[compared == TRUE & match == FALSE, .(dgp, N, ntime, label, sim_id, used_saved, path_pkg, code_saved, est_saved, est_pkg, se_saved, se_pkg)]) }
cat("largest relative |est| and |se| difference (numeric comparisons):",
    R[numeric_compared == TRUE, max(abs(est_pkg - est_saved) / abs(est_saved))],
    R[numeric_compared == TRUE, max(abs(se_pkg - se_saved) / abs(se_saved))], "\n")
cat("description: default control (check_coherence = TRUE) gives a different reported fit in", R[, sum(default_differs)], "of",
    R[pkg_status %in% "fitted" & default_status %in% "fitted", .N], "fitted datasets\n")
cat("study rows read unchanged (md5 before = after):", identical(md5_before, md5_after), "\n")
GATE3 <- R[compared == TRUE, all(match)] && R[, sum(compared)] >= 20L && R[generated %in% TRUE, all(data_identical %in% TRUE)] &&
  R[, all(is.na(error))] && identical(md5_before, md5_after)
cat("GATE check 3:", if (GATE3) "PASS" else "FAIL", "\n")
sink()
