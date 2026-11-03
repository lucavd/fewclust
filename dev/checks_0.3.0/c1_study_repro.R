# c1_study_repro.R -- fewclust 0.3.0, checks 1 and 2 (pre-registered in dev/aggiornamento_0.3.0.md).
# Regenerates study datasets from their seeds with the code of branch morris-v4 (read only) and compares
#   check 1: flip_test() with the saved p-value of the paper's sign-flip test (p0_guard of results_v4/S7d2/rows,
#            the version with the non-computability rules of Section 2.4), same seed 20261022 + i, same signs;
#   check 2: auto_cluster_fit() with the saved GEE AR(1) + Fay-Graubard estimate and standard error of the
#            same replicate (ar1_est, ar1_se_fg in the same rows; geeglm scale.fix = TRUE, fewclust 0.2.0 FG).
# Usage (from the package root): Rscript --vanilla dev/checks_0.3.0/c1_study_repro.R <dev_lib> <out_dir>
# Writes only to <out_dir>.
args <- commandArgs(TRUE)
DEV <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PROJ <- "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
.libPaths(c(DEV, PROJ))
invisible(loadNamespace("fewclust"))              # 0.3.0 from the dev library, before the study code loads 0.2.0
stopifnot(getNamespaceVersion("fewclust") == "0.3.0")
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))  # generators, seeds, tasks, flip_prep/flip_signs/flip_p_b
# the study code calls library(fewclust), which under R 4.6 replaces the loaded 0.3.0 with the 0.2.0 of ~/R/s7d_lib:
# detach it and load the development version again (the study functions used here do not call fewclust)
if ("package:fewclust" %in% search()) detach("package:fewclust", unload = TRUE)
if ("fewclust" %in% loadedNamespaces()) unloadNamespace("fewclust")
invisible(loadNamespace("fewclust", lib.loc = DEV))
stopifnot(getNamespaceVersion("fewclust") == "0.3.0", dirname(getNamespaceInfo("fewclust", "path")) == DEV)
suppressMessages(library(data.table)); setDTthreads(1L)

TA <- s7d2_tasks()
MECH <- c("GenBinaryY", "GLMM-native", "GLMM-AR1", "Copula-AR1")
stopifnot(all(MECH %in% TA$dgp))
pick <- function(dg, sc, n, nt, lb, ids, set) {
  r <- which(TA$dgp == dg & TA$scen == sc & TA$N == n & TA$ntime == nt & TA$label == lb & TA$chunk == 1L)
  stopifnot(length(r) == 1L)
  data.table(row = r, sim_id = ids, set = set)
}
# A1: scenario M1, four mechanisms, N in {10, 25}, T in {2, 10}, H0 and H1, replicates 1-3
A1 <- rbindlist(lapply(MECH, function(dg) rbindlist(lapply(c(10L, 25L), function(n) rbindlist(lapply(c(2L, 10L), function(nt)
  rbindlist(lapply(c("manual", "mmlb"), function(lb) pick(dg, "M1", n, nt, lb, 1:3, "A1")))))))))
# A2: covariates in the mean model (M2: age; M5: dose), H0, replicates 1-3
A2 <- rbind(pick("GenBinaryY", "M2", 25L, 3L, "manual", 1:3, "A2"), pick("GLMM-native", "M5", 10L, 3L, "manual", 1:3, "A2"))
# A3: datasets where the paper's test was not computable: M1, N = 10, T = 2, four mechanisms, H0 and H1,
#     the first three replicates after 3 with a saved reason (why0_guard)
A3 <- rbindlist(lapply(MECH, function(dg) rbindlist(lapply(c("manual", "mmlb"), function(lb) {
  r <- pick(dg, "M1", 10L, 2L, lb, NA_integer_, "A3")$row
  s <- readRDS(TA$file[r])
  ids <- head(s[sim_id > 3L & !is.na(why0_guard), sim_id], 3L)
  if (length(ids)) data.table(row = r, sim_id = ids, set = "A3") else NULL
}))))
POOL <- rbind(A1, A2, A3)
cat("datasets in the pool:", nrow(POOL), "(A1", nrow(A1), "A2", nrow(A2), "A3", nrow(A3), ")\n")

files <- unique(TA$file[POOL$row])
md5_before <- tools::md5sum(files)

one <- function(k) {
  t <- TA[POOL$row[k]]; i <- POOL$sim_id[k]
  gen <- s7d2_gen(t)
  saved <- readRDS(t$file)[sim_id == i]
  stopifnot(nrow(saved) == 1L)
  deg <- design_degenerate(t$model, t$N, t$ntime, i)   # before the seeded generation, as in one_rep_s7d2
  set_rep_seed(t$model, i)
  d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) e)
  if (inherits(d, "error")) {   # as in one_rep_s7d2: no dataset, nothing to compare but the saved flag
    return(data.table(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label, sim_id = i, set = POOL$set[k],
                      design_degenerate = deg, generated = FALSE, generated_saved = saved$generated,
                      data_identical = isFALSE(saved$generated), gen_error = substr(conditionMessage(d), 1, 120)))
  }
  d <- d[order(d$id, d$time), ]
  form <- FORM_OF(t$scen)
  term <- tail(colnames(model.matrix(form, d)), 1L)
  stopifnot(term == "time:trt")
  # check 1: the package's test and, for diagnosis, the study code run live
  ft <- fewclust::flip_test(form, d, group = "id", time = "time", term = term, seed = FLIP_SEED(i))
  set.seed(FLIP_SEED(i)); P <- flip_prep(d, t$scen); S <- flip_signs(P$N); live <- flip_p_b(P, S, 0)
  # check 2: the automatic strategy of the package
  fit <- tryCatch(fewclust::auto_cluster_fit(form, d, group = "id", time = "time"), error = function(e) e)
  if (inherits(fit, "error")) fit <- structure(list(decision = list(rule = paste("error:", conditionMessage(fit)))),
                                               class = "fewclust_refusal")
  refused <- inherits(fit, "fewclust_refusal")
  j <- if (!refused) match(term, fit$coefficients$term) else NA_integer_
  data.table(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label, sim_id = i, set = POOL$set[k],
             design_degenerate = deg, generated = TRUE, generated_saved = saved$generated,
             data_identical = nrow(d) == saved$n_obs && sum(d$binY) == saved$sum_y,
             p_saved = saved$p0_guard, why_saved = saved$why0_guard,
             p_pkg = ft$p.value, why_pkg = ft$reason, n_signs = ft$n_signs, p_live = live$p,
             gee_saved_ok = isTRUE(saved$ar1_gee_ok) && is.finite(saved$ar1_se_fg),
             est_saved = saved$ar1_est, se_saved = saved$ar1_se_fg,
             pkg_status = if (refused) paste("refused:", fit$decision$rule) else "fitted",
             pkg_strategy = if (refused) NA_character_ else fit$strategy,
             est_pkg = if (refused) NA_real_ else fit$coefficients$estimate[j],
             se_pkg = if (refused) NA_real_ else fit$coefficients$std.error[j])
}
one_safe <- function(k) tryCatch(one(k), error = function(e) {
  t <- TA[POOL$row[k]]
  data.table(dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label, sim_id = POOL$sim_id[k],
             set = POOL$set[k], error = conditionMessage(e))
})
t0 <- proc.time()[[3]]
R <- rbindlist(parallel::mclapply(seq_len(nrow(POOL)), one_safe, mc.cores = 4L, mc.preschedule = FALSE), fill = TRUE)
if (!"error" %in% names(R)) R[, error := NA_character_]
cat("seconds:", round(proc.time()[[3]] - t0), "\n")
stopifnot(nrow(R) == nrow(POOL))

cat("datasets with an error in the check:", R[!is.na(error), .N], "\n"); print(R[!is.na(error), .(dgp, N, ntime, label, sim_id, set, error)])
R[, flip_match := (is.na(p_saved) & is.na(p_pkg) & why_saved == why_pkg) |
                  (!is.na(p_saved) & !is.na(p_pkg) & abs(p_saved - p_pkg) <= 1e-12)]
# a dataset that the study code cannot generate (as in the saved run) has nothing to compare
R[generated == FALSE, flip_match := generated_saved %in% FALSE]
R[is.na(flip_match), flip_match := FALSE]
R[, gee_compared := !is.na(pkg_status) & pkg_status == "fitted" & gee_saved_ok %in% TRUE]
R[, gee_match := ifelse(gee_compared, pkg_strategy == "gee_ar1_fg" &
                          abs(est_pkg - est_saved) <= 1e-8 * pmax(1, abs(est_saved)) &
                          abs(se_pkg - se_saved) <= 1e-8 * pmax(1, abs(se_saved)), NA)]
md5_after <- tools::md5sum(files)

fwrite(R, file.path(OUT, "c1_study_repro.csv"))
sink(file.path(OUT, "c1_summary.txt"), split = TRUE)
cat("fewclust", as.character(getNamespaceVersion("fewclust")), "\n")
cat("datasets:", nrow(R), "; with an error in the check:", R[!is.na(error), .N],
    "; data identical to the saved run (n_obs, sum_y):", R[, sum(data_identical, na.rm = TRUE)], "\n")
if (R[!is.na(error), .N]) print(R[!is.na(error), .(dgp, N, ntime, label, sim_id, set, error)])
cat("check 1 (sign-flip p-value equal to the saved one, |diff| <= 1e-12, or both not computable with the same reason):",
    R[, sum(flip_match)], "of", nrow(R), "\n")
cat("  with enumeration of all signs:", R[n_signs < 2000, .N], "; with 1,999 random signs + identity:", R[n_signs == 2000, .N],
    "; not computable (saved):", R[generated == TRUE & is.na(p_saved), .N],
    "; not generated, here and in the saved run:", R[generated == FALSE & generated_saved %in% FALSE, .N], "\n")
print(R[generated == TRUE & is.na(p_saved), .N, by = .(why_saved, why_pkg)])
cat("  degenerate designs (excluded from every denominator of the paper):", R[design_degenerate == TRUE, .N], "\n")
cat("check 2 (auto_cluster_fit = saved GEE AR(1) + Fay-Graubard, relative tolerance 1e-8):",
    R[gee_compared == TRUE, sum(gee_match)], "of", R[, sum(gee_compared)], "compared\n")
cat("  not compared: refused by the package", R[!is.na(pkg_status) & pkg_status != "fitted", .N],
    "; saved GEE or FG not usable", R[pkg_status %in% "fitted" & !(gee_saved_ok %in% TRUE), .N], "\n")
print(R[gee_compared == FALSE, .N, by = .(pkg_status, pkg_strategy, gee_saved_ok)])
cat("largest |p_pkg - p_saved|:", R[!is.na(p_saved) & !is.na(p_pkg), max(abs(p_pkg - p_saved))], "\n")
cat("largest relative |est| and |se| difference:", R[gee_compared == TRUE, max(abs(est_pkg - est_saved) / pmax(1, abs(est_saved)))],
    R[gee_compared == TRUE, max(abs(se_pkg - se_saved) / pmax(1, abs(se_saved)))], "\n")
cat("study rows read unchanged (md5 before = after):", identical(md5_before, md5_after), "\n")
GATE1 <- R[, all(flip_match)] && R[, all(data_identical %in% TRUE)] && nrow(R) >= 20L
GATE2 <- R[gee_compared == TRUE, all(gee_match)] && R[, sum(gee_compared)] >= 20L
cat("GATE check 1:", if (GATE1) "PASS" else "FAIL", "\nGATE check 2:", if (GATE2) "PASS" else "FAIL", "\n")
sink()
