# c4_paper_condition.R -- fewclust 0.3.1, check 10: the condition pre-registered in the paper for citing fewclust as
# the implementation of its rule (branch morris-v4, R/s7f/s7f_fewclust_condition.R, decisions D-43 and D-65), run with
# the same datasets and the same code, against the development library, with the condition adapted to the rule of the
# corrected study (D-82): the chosen strategy is the one of tables/V5_S7f_table4.csv (column chosen), i.e. "mmm".
# The script on branch morris-v4 still looks for "ar1" and "fg" and must be updated there. Writes only to <out_dir>.
# Usage: Rscript --vanilla dev/checks_0.3.1/c4_paper_condition.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages({ library(fewclust); library(data.table) })
stopifnot(as.character(packageVersion("fewclust")) == "0.3.1")
T4 <- read.csv(path.expand("~/MMM/tables/V5_S7f_table4.csv")); rule <- c(mmm_cc = "mmm", gee_ar1_fg = "gee_ar1_fg")[[T4$strategy[T4$chosen]]]
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
R[, paper_rule := rule]; R[, condition_met := chosen == rule]
fwrite(R, file.path(OUT, "c4_paper_condition.csv"))
sink(file.path(OUT, "c4_summary.txt"), split = TRUE); print(R)
cat("GATE check 10:", if (all(R$condition_met)) "PASS" else "FAIL", "\n"); sink()
