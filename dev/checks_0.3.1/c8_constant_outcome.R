# c8_constant_outcome.R -- fewclust 0.3.1, check 8 (pre-registered in dev/aggiornamento_0.3.1.md): a dataset whose
# outcome is constant (all 0 or all 1) gives no usable result for any strategy, as decided for the study (D-82 on branch
# morris-v4). 16 constructed datasets: simulate_clustered_binary() with N in {10, 25}, T in {2, 4}, with and without a
# covariate, outcome set to all 0 and to all 1. On each: cluster_strategies() refuses with rule "D2: constant outcome";
# auto_cluster_fit() returns a fewclust_refusal with that rule for "auto" and for each of the 19 accepted strategy names;
# flip_test() returns p = NA, computable = FALSE, a reason starting with "constant outcome". Counter-check: with one
# outcome changed the reason "constant outcome" does not appear. Usage: Rscript --vanilla dev/checks_0.3.1/c8_constant_outcome.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(fewclust)); stopifnot(as.character(packageVersion("fewclust")) == "0.3.1")
S <- c("auto", fewclust:::fc_strategy_names()); stopifnot(length(S) == 20L)
G <- expand.grid(N = c(10L, 25L), T = c(2L, 4L), cov = c(FALSE, TRUE), val = 0:1)
rows <- list(); ctrl <- list()
for (k in seq_len(nrow(G))) {
  g <- G[k, ]
  d <- simulate_clustered_binary(n_groups = g$N, n_occasions = g$T, seed = k)
  d$y <- g$val
  f <- if (g$cov) y ~ time * arm + x else y ~ time * arm
  pl <- cluster_strategies(f, d, group = id, time = time)
  fits <- vapply(S, function(s) {
    r <- tryCatch(auto_cluster_fit(f, d, group = id, time = time, strategy = s), error = function(e) e)
    inherits(r, "fewclust_refusal") && identical(r$decision$rule, "D2: constant outcome") }, logical(1))
  ft <- flip_test(f, d, group = id, time = time, seed = k)
  rows[[k]] <- data.frame(N = g$N, T = g$T, covariate = g$cov, outcome = g$val,
    plan_refused_D2 = pl$refused && pl$rule == "D2: constant outcome", fits_refused_D2 = sum(fits), n_strategies = length(S),
    flip_na = is.na(ft$p.value), flip_not_computable = !ft$computable,
    flip_reason = ft$reason, flip_reason_ok = isTRUE(startsWith(ft$reason, "constant outcome")))
  d2 <- d; d2$y[1] <- 1L - g$val
  f2 <- flip_test(f, d2, group = id, time = time, seed = k)
  ctrl[[k]] <- data.frame(N = g$N, T = g$T, covariate = g$cov, outcome = g$val, p = f2$p.value, reason = f2$reason,
                          no_constant_reason = !isTRUE(startsWith(f2$reason, "constant outcome")))
}
R <- do.call(rbind, rows); C <- do.call(rbind, ctrl)
R$ok <- R$plan_refused_D2 & R$fits_refused_D2 == R$n_strategies & R$flip_na & R$flip_not_computable & R$flip_reason_ok
write.csv(R, file.path(OUT, "c8_constant_outcome.csv"), row.names = FALSE)
write.csv(C, file.path(OUT, "c8_counter_check.csv"), row.names = FALSE)
sink(file.path(OUT, "c8_summary.txt"), split = TRUE)
print(R[, c("N", "T", "covariate", "outcome", "plan_refused_D2", "fits_refused_D2", "flip_na", "flip_reason_ok", "ok")], row.names = FALSE)
cat("flip_test reason:", unique(R$flip_reason), "\n")
cat("counter-check (one outcome changed): reason 'constant outcome' absent in", sum(C$no_constant_reason), "of", nrow(C), "\n")
cat("GATE check 8:", if (all(R$ok) && all(C$no_constant_reason) && nrow(R) == 16L) "PASS" else "FAIL",
    "(", sum(R$ok), "of", nrow(R), "datasets;", sum(R$fits_refused_D2), "of", sum(R$n_strategies), "fits refused with D2 )\n")
sink()
