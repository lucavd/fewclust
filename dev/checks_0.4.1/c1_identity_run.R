# c1_identity_run.R -- fewclust 0.4.1, gates A and B on the 60 datasets of check 1 of 0.4.0 (dev/aggiornamento_0.4.1.md).
# Run once with the 0.4.0 reference library (it writes the datasets: simulate_clustered_binary(), seeds 1-60, as in
# dev/checks_0.4.0/c1_identity_run.R) and once with the 0.4.1 library (it reads them). For each dataset: the 14 calls of
# check 1 of 0.4.0 (auto, mmm, GEE 3 working correlations x 4 corrections; default control), plus auto and mmm with
# fewclust_control(check_coherence = FALSE); for each call the table, or the error, and the record of the convergence
# check (fit$convergence, or the field `convergence` of the error in 0.4.1); flip_test() with seed 100 + s; gee_vcov().
# Usage: Rscript --vanilla c1_identity_run.R <lib> <datasets.rds> <out.rds>
args <- commandArgs(TRUE); LIB <- args[1]; DS <- args[2]; OUT <- args[3]
.libPaths(c(LIB, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(fewclust, lib.loc = LIB); library(geepack) })
ver <- as.character(packageVersion("fewclust", lib.loc = LIB)); stopifnot(ver %in% c("0.4.0", "0.4.1"))
if (ver == "0.4.0") {
  grid <- expand.grid(N = c(10L, 15L, 25L, 65L), T = c(2L, 3L, 4L, 8L), rep = 1:4)
  grid <- grid[order(grid$rep, grid$N, grid$T), ][1:60, ]
  D <- lapply(seq_len(nrow(grid)), function(s) {
    d <- as.data.frame(simulate_clustered_binary(n_groups = grid$N[s], n_occasions = grid$T[s],
                                                 slope_sd = if (s %% 3 == 0) 0.5 else 0, seed = s))
    list(data = d[order(d$id, d$time), ], formula = if (s %% 2 == 0) y ~ time * arm + x else y ~ time * arm)
  })
  saveRDS(D, DS)
} else D <- readRDS(DS)
cases <- rbind(data.frame(method = c("auto", "mmm"), corstr = "ar1", correction = "fg", coherence = TRUE),
               expand.grid(method = "gee", corstr = c("ar1", "exchangeable", "independence"),
                           correction = c("fg", "md", "kc", "none"), coherence = TRUE, stringsAsFactors = FALSE),
               data.frame(method = c("auto", "mmm"), corstr = "ar1", correction = "fg", coherence = FALSE))
run_case <- function(d, f, cs) {
  r <- tryCatch(suppressWarnings(suppressMessages(
    fit_fewclust(d, f, id = id, time = time, method = cs$method, corstr = cs$corstr, correction = cs$correction,
                 control = fewclust_control(check_coherence = cs$coherence)))), fewclust_error = function(e) e)
  if (inherits(r, "error")) return(list(refused = TRUE, message = conditionMessage(r), convergence = r$convergence))
  list(refused = FALSE, method = r$method, table = as.data.frame(tidy(r)), fallback = r$fallback,
       convergence = r$convergence, code = if (r$method == "mmm") r$engine$code else NA_integer_)
}
res <- lapply(seq_along(D), function(s) {
  d <- D[[s]]$data; f <- D[[s]]$formula
  fits <- lapply(seq_len(nrow(cases)), function(k) run_case(d, f, cases[k, ]))
  fl <- suppressWarnings(flip_test(d, f, id = id, time = time, seed = 100 + s))
  vc <- list()
  for (cs in c("ar1", "exchangeable", "independence")) for (sf in c(FALSE, TRUE)) {
    m <- tryCatch(geeglm(f, id = id, data = d, family = binomial, corstr = cs, scale.fix = sf), error = function(e) NULL)
    vc[[paste(cs, sf)]] <- if (is.null(m)) "geeglm error" else
      lapply(c("unc", "fg", "md", "kc"), function(ty) tryCatch(gee_vcov(m, ty), error = conditionMessage))
  }
  list(fits = fits, flip = c(p = fl$p.value, reason = fl$reason), vcov = vc)
})
saveRDS(list(version = ver, cases = cases, res = res), OUT)
cat("fewclust", ver, ": datasets", length(res), "; fits", length(res) * nrow(cases), "\n")
