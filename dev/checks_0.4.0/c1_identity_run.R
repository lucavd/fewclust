# c1_identity_run.R -- fewclust 0.4.0, checks 1 and 2 (dev/aggiornamento_0.4.0.md). Run once with the 0.3.1
# reference library (it also writes the 60 datasets) and once with the development library (it reads them).
# For each dataset: every strategy of 0.3.1 mapped to the 0.4.0 call, the coefficient table or the refusal,
# flip_test() with seed 100 + s, and gee_vcov() on geeglm fits (3 working correlations x scale.fix x 4 types).
# Usage: Rscript --vanilla c1_identity_run.R <lib> <datasets.rds> <out.rds>
args <- commandArgs(TRUE); LIB <- args[1]; DS <- args[2]; OUT <- args[3]
.libPaths(c(LIB, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(fewclust); library(geepack) })
ver <- as.character(packageVersion("fewclust")); old <- ver == "0.3.1"
stopifnot(ver %in% c("0.3.1", "0.4.0"))
if (old) {
  grid <- expand.grid(N = c(10L, 15L, 25L, 65L), T = c(2L, 3L, 4L, 8L), rep = 1:4)
  grid <- grid[order(grid$rep, grid$N, grid$T), ][1:60, ]
  D <- lapply(seq_len(nrow(grid)), function(s) {
    d <- as.data.frame(simulate_clustered_binary(n_groups = grid$N[s], n_occasions = grid$T[s],
                                                 slope_sd = if (s %% 3 == 0) 0.5 else 0, seed = s))
    list(data = d[order(d$id, d$time), ], formula = if (s %% 2 == 0) y ~ time * arm + x else y ~ time * arm)
  })
  saveRDS(D, DS)
} else D <- readRDS(DS)
tag <- c(ar1 = "ar1", exchangeable = "exch", independence = "ind")
corr <- c(fg = "fg", md = "md", kc = "kc", none = "unc")
cases <- rbind(data.frame(method = c("auto", "mmm"), corstr = "ar1", correction = "fg"),
               expand.grid(method = "gee", corstr = names(tag), correction = names(corr), stringsAsFactors = FALSE))
run_case <- function(d, f, cs) {
  if (old) {
    strat <- if (cs$method == "gee") paste0("gee_", tag[[cs$corstr]], "_", corr[[cs$correction]]) else cs$method
    r <- auto_cluster_fit(f, d, group = id, time = time, strategy = strat)
    if (inherits(r, "fewclust_refusal")) return(list(refused = TRUE, rule = r$decision$rule))
    list(refused = FALSE, method = if (startsWith(r$strategy, "gee")) "gee" else "mmm",
         table = as.data.frame(r$coefficients)[, c("term", "estimate", "std.error", "statistic", "df", "p.value")])
  } else {
    r <- tryCatch(suppressWarnings(suppressMessages(
      fit_fewclust(d, f, id = id, time = time, method = cs$method, corstr = cs$corstr,
                   correction = cs$correction))), fewclust_error = function(e) e)
    if (inherits(r, "error")) return(list(refused = TRUE, rule = conditionMessage(r)))
    list(refused = FALSE, method = r$method, table = as.data.frame(tidy(r)))
  }
}
res <- lapply(seq_along(D), function(s) {
  d <- D[[s]]$data; f <- D[[s]]$formula
  fits <- lapply(seq_len(nrow(cases)), function(k) run_case(d, f, cases[k, ]))
  fl <- if (old) suppressWarnings(flip_test(f, d, group = id, time = time, seed = 100 + s)) else
    suppressWarnings(flip_test(d, f, id = id, time = time, seed = 100 + s))
  vc <- list()
  for (cs in names(tag)) for (sf in c(FALSE, TRUE)) {
    m <- tryCatch(geeglm(f, id = id, data = d, family = binomial, corstr = cs, scale.fix = sf), error = function(e) NULL)
    vc[[paste(cs, sf)]] <- if (is.null(m)) "geeglm error" else
      lapply(c("unc", "fg", "md", "kc"), function(ty) tryCatch(gee_vcov(m, ty), error = conditionMessage))
  }
  list(fits = fits, flip = c(p = fl$p.value, reason = fl$reason), vcov = vc)
})
saveRDS(list(version = ver, cases = cases, res = res), OUT)
cat("fewclust", ver, ": datasets", length(res), "; fits", length(res) * nrow(cases), "\n")
