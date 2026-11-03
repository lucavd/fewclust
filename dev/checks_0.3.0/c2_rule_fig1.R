# c2_rule_fig1.R -- fewclust 0.3.0, check 3: the rule of the package gives the choice of Figure 1 of the paper
# (manuscript v5, branch morris-v4, py/s7f/fig1_rule_v5.json) for N in {10, 25, 145} and T in {2, 3, 10}.
# Usage: Rscript --vanilla dev/checks_0.3.0/c2_rule_fig1.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(fewclust)); stopifnot(as.character(packageVersion("fewclust")) == "0.3.0")
fig1 <- readLines(path.expand("~/MMM/py/s7f/fig1_rule_v5.json"))
cat("Figure 1 boxes (fig1_rule_v5.json):\n"); cat(fig1[grepl("GEE|marginalized|occasions|sign-flip|closer|clusters", fig1)], sep = "\n")
# expected choice, transcribed from the four boxes of Figure 1:
#   fg  : GEE with AR(1) working correlation, Fay-Graubard correction, t reference (N - p df); any N and T
#   mmm : alternative, marginalized multilevel model, t reference, convergence check; closer to 5% with 2 occasions and 15+ clusters
#   end : 2-3 occasions, <= 25 clusters: both answer in fewer than half the datasets; sign-flip score test as a complement
G <- expand.grid(N = c(10L, 25L, 145L), T = c(2L, 3L, 10L))
rows <- lapply(seq_len(nrow(G)), function(k) {
  N <- G$N[k]; T <- G$T[k]
  d <- simulate_clustered_binary(n_groups = N, n_occasions = T, seed = 1)
  pl <- cluster_strategies(y ~ time * arm, d, group = id, time = time)
  data.frame(N = N, T = T, refused = pl$refused, strategy = pl$strategy, reference = pl$reference,
             alternative = pl$alternative, complement = pl$complement,
             mmm_closer_note = any(grepl("slightly closer to 5%", pl$cautions)),
             exp_strategy = "gee_ar1_fg", exp_reference = "t", exp_alternative = "mmm",
             exp_complement = if (T <= 3 && N <= 25) "flip_test" else NA_character_,
             exp_mmm_closer_note = T == 2 && N >= 15)
})
R <- do.call(rbind, rows)
R$match <- !R$refused & R$strategy == R$exp_strategy & R$reference == R$exp_reference &
  R$alternative == R$exp_alternative & (R$complement %in% R$exp_complement) & R$mmm_closer_note == R$exp_mmm_closer_note
write.csv(R, file.path(OUT, "c2_rule_fig1.csv"), row.names = FALSE)
sink(file.path(OUT, "c2_summary.txt"), split = TRUE)
print(R[, c("N", "T", "strategy", "reference", "alternative", "complement", "mmm_closer_note", "match")], row.names = FALSE)
cat("GATE check 3:", if (all(R$match)) "PASS" else "FAIL", "(", sum(R$match), "of", nrow(R), "cells as in Figure 1 )\n")
sink()
