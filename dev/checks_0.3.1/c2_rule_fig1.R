# c2_rule_fig1.R -- fewclust 0.3.1, check 4: the rule of the package on the 3 x 3 grid N in {10, 25, 145} x T in {2, 3, 10}
# gives the choice of the corrected study (decision D-82 of branch morris-v4): the strategy chosen in tables/V5_S7f_table4.csv
# (column chosen), the GEE AR(1) + Fay-Graubard alternative, the sign-flip complement with T <= 3 and N <= 25, and a region
# note exactly where tables/V5_S7f_h2h_boot.csv (Figure 1a, q = mad_gee_minus_mmm, regions map|...) has class GEE or MMM.
# Expected values are read from the two files, not transcribed. Usage: Rscript --vanilla dev/checks_0.3.1/c2_rule_fig1.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(fewclust)); stopifnot(as.character(packageVersion("fewclust")) == "0.3.1")
T4 <- read.csv(path.expand("~/MMM/tables/V5_S7f_table4.csv"))
H <- read.csv(path.expand("~/MMM/tables/V5_S7f_h2h_boot.csv"))
chosen <- T4$strategy[T4$chosen]; stopifnot(length(chosen) == 1L)
exp_strategy <- c(mmm_cc = "mmm", gee_ar1_fg = "gee_ar1_fg")[[chosen]]
exp_alternative <- setdiff(c("mmm", "gee_ar1_fg"), exp_strategy)
cls <- H[grepl("^map\\|", H$region) & H$q == "mad_gee_minus_mmm", c("region", "class")]
cat("chosen in table4:", chosen, "; Figure 1a classes:\n"); print(cls, row.names = FALSE)
band <- function(N, T) sprintf("map|%s|%s", if (N < 15) "10" else if (N <= 25) "15-25" else ">25",
                               if (T <= 2) "2" else if (T <= 3) "3" else "4+")
G <- expand.grid(N = c(10L, 25L, 145L), T = c(2L, 3L, 10L))
rows <- lapply(seq_len(nrow(G)), function(k) {
  N <- G$N[k]; T <- G$T[k]
  d <- simulate_clustered_binary(n_groups = N, n_occasions = T, seed = 1)
  pl <- cluster_strategies(y ~ time * arm, d, group = id, time = time)
  cl <- cls$class[cls$region == band(N, T)]; stopifnot(length(cl) == 1L)
  data.frame(N = N, T = T, region = band(N, T), class = cl, refused = pl$refused, strategy = pl$strategy,
             reference = pl$reference, alternative = pl$alternative, complement = pl$complement,
             region_note = any(grepl("closer to 5% than", pl$cautions)),
             old_note = any(grepl("slightly closer to 5%", pl$cautions)),
             exp_strategy = exp_strategy, exp_alternative = exp_alternative,
             exp_complement = if (T <= 3 && N <= 25) "flip_test" else NA_character_,
             exp_region_note = cl %in% c("GEE", "MMM"))
})
R <- do.call(rbind, rows)
R$match <- !R$refused & R$strategy == R$exp_strategy & R$reference == "t" & R$alternative == R$exp_alternative &
  (R$complement %in% R$exp_complement) & R$region_note == R$exp_region_note & !R$old_note
write.csv(R, file.path(OUT, "c2_rule_fig1.csv"), row.names = FALSE)
sink(file.path(OUT, "c2_summary.txt"), split = TRUE)
print(R[, c("N", "T", "class", "strategy", "reference", "alternative", "complement", "region_note", "match")], row.names = FALSE)
cat("GATE check 4:", if (all(R$match)) "PASS" else "FAIL", "(", sum(R$match), "of", nrow(R), "cells as in the corrected study )\n")
sink()
