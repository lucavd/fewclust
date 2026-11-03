# c4_paper_condition.R -- fewclust 0.3.0, check 7: the condition pre-registered in the paper for citing fewclust as
# the implementation of its rule (branch morris-v4, R/s7f/s7f_fewclust_condition.R, decisions D-43 and D-65), run with
# the same datasets and the same test, against the development library. Writes only to <out_dir>.
# Usage: Rscript --vanilla dev/checks_0.3.0/c4_paper_condition.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(fewclust); library(data.table) })
stopifnot(as.character(packageVersion("fewclust")) == "0.3.0")
data(fc_trial, package = "fewclust")
sim <- simulate_clustered_binary(n_groups = 25, n_occasions = 4, seed = 20261027)
# from here to the condition: the code of R/s7f/s7f_fewclust_condition.R, unchanged
pick <- function(nm, cand) cand[cand %in% nm][1]
run <- function(d, label) {
  nm <- names(d); y <- pick(nm, c("y", "outcome", "resp")); g <- pick(nm, c("id", "cluster", "group", "subject"))
  tt <- pick(nm, c("time", "occasion", "visit")); tr <- pick(nm, c("trt", "treatment", "arm", "x"))
  f <- as.formula(sprintf("%s ~ %s * %s", y, tr, tt))
  fit <- tryCatch(auto_cluster_fit(f, data = d, group = g, time = tt), error = function(e) e)
  chosen <- if (inherits(fit, "error")) paste("error:", conditionMessage(fit)) else {
    s <- fit$strategy; if (is.null(s)) s <- attr(fit, "strategy"); if (is.null(s)) s <- class(fit)[1]; paste(unlist(s), collapse = "/") }
  data.table(dataset = label, formula = deparse(f), group = g, time = tt, chosen = chosen)
}
R <- rbind(run(as.data.frame(fc_trial), "fc_trial"), run(as.data.frame(sim), "simulated N=25 T=4 seed 20261027"))
R[, paper_rule := "gee_ar1_fg_t"]; R[, condition_met := grepl("ar1", chosen, ignore.case = TRUE) & grepl("fg|fay", chosen, ignore.case = TRUE)]
fwrite(R, file.path(OUT, "c4_paper_condition.csv"))
sink(file.path(OUT, "c4_summary.txt"), split = TRUE); print(R)
cat("GATE check 7:", if (all(R$condition_met)) "PASS" else "FAIL", "\n"); sink()
