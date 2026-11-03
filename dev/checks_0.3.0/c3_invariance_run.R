# c3_invariance_run.R -- fewclust 0.3.0, check 6 (part 1): computes, with the fewclust found first in <lib>,
# the quantities the paper took from fewclust (gee_vcov: uncorrected, FG, MD, KC on geeglm fits with the three
# working correlations, scale.fix TRUE and FALSE) and the marginalized-model engine (fit_mmm with an explicit
# control), and saves them to <out_rds>. Run once with 0.2.0 (~/R/s7d_lib) and once with 0.3.0; c3_invariance_compare.R
# compares the two files with identical().
# Usage: Rscript --vanilla dev/checks_0.3.0/c3_invariance_run.R <lib> <out_rds>
args <- commandArgs(TRUE); LIB <- normalizePath(args[1]); OUT <- args[2]
.libPaths(c(LIB, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(fewclust); library(geepack) })
ver <- as.character(packageVersion("fewclust"))
sim <- function(s) {
  N <- c(10L, 25L)[1 + s %% 2]; T <- c(2L, 4L, 6L)[1 + s %% 3]
  d <- simulate_clustered_binary(n_groups = N, n_occasions = T, seed = s)
  d[order(d$id, d$time), ]
}
V <- list()
for (s in 1:100) {
  d <- sim(s)
  for (cs in c("exchangeable", "independence", "ar1")) for (sf in c(FALSE, TRUE)) {
    m <- tryCatch(geeglm(y ~ time * arm, id = id, data = d, family = binomial, corstr = cs, scale.fix = sf),
                  error = function(e) NULL)
    key <- sprintf("s%03d_%s_%s", s, cs, sf)
    V[[key]] <- if (is.null(m)) "geeglm error" else lapply(c(unc = "unc", fg = "fg", md = "md", kc = "kc"), function(ty)
      tryCatch(gee_vcov(m, ty), error = function(e) conditionMessage(e)))
  }
}
M <- list()
fx_ctrl <- fewclust_control(corstr = "exchangeable", check_coherence = TRUE)
for (s in 1:40) {
  d <- sim(1000 + s)
  p <- fewclust:::fc_parse(y ~ time * arm, d, quote(id), quote(time))
  r <- fewclust:::fit_mmm(p, fx_ctrl, "t", 0.05, nlevels(p$id) - ncol(p$X))
  M[[sprintf("m%03d", s)]] <- if (inherits(r, "fc_failure")) r$reason else
    list(coef = r$coefficients, vcov = r$vcov, theta = r$model$theta, code = r$model$code,
         boundary = r$model$check$boundary, latent_sd = r$model$check$latent_sd)
}
saveRDS(list(version = ver, V = V, M = M), OUT)
cat("fewclust", ver, ": gee_vcov sets", length(V), "; marginalized-model fits", length(M), "\n")
