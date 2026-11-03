# c8_constant_outcome.R -- fewclust 0.3.0, observation (not pre-registered): how many datasets of the study with a
# constant outcome (all 0 or all 1) count as usable for GEE AR(1) + Fay-Graubard under the usability filter of the paper
# (Section 2.5: finite estimate, SE and p; |estimate| < 10; 0 < SE < 10), and how much they add to its usable share
# (Table 4: mean over design points of the share of datasets with an estimable design). Reads the saved rows of
# results_v4/S7d2/rows (branch morris-v4, read only). Writes only to <out_dir>.
# Usage: Rscript --vanilla dev/checks_0.3.0/c8_constant_outcome.R <out_dir>
args <- commandArgs(TRUE); OUT <- args[1]
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(1L)
files <- list.files(path.expand("~/MMM/results_v4/S7d2/rows"), "\\.rds$", full.names = TRUE)
stopifnot(length(files) == 47520L)
md5_before <- tools::md5sum(files[c(1, 23760, 47520)])
cols <- c("dgp", "scen", "N", "ntime", "label", "sim_id", "design_degenerate", "n_obs", "sum_y", "ar1_gee_ok", "ar1_est", "ar1_se_fg")
X <- rbindlist(parallel::mclapply(files, function(f) readRDS(f)[, ..cols], mc.cores = 4L))
X <- X[design_degenerate == FALSE]
X[, usable := ar1_gee_ok %in% TRUE & is.finite(ar1_est) & is.finite(ar1_se_fg) & abs(ar1_est) < 10 & ar1_se_fg > 0 & ar1_se_fg < 10]
X[, constant := !is.na(sum_y) & (sum_y == 0 | sum_y == n_obs)]
D <- X[, .(n = .N, pu = 100 * mean(usable), pu_noconst = 100 * mean(usable & !constant),
           const = 100 * mean(constant), const_usable = sum(constant & usable)), by = .(dgp, scen, N, ntime, label)]
H0 <- D[label == "manual"]
stopifnot(nrow(H0) == 2376L)
small <- H0[ntime <= 3 & N <= 25]
T4 <- read.csv(path.expand("~/MMM/tables/V5_S7f_table4.csv"))
fwrite(D, file.path(OUT, "c8_constant_outcome_dp.csv"))
sink(file.path(OUT, "c8_summary.txt"), split = TRUE)
cat("replicates with an estimable design:", nrow(X), "; with a constant outcome:", X[, sum(constant)],
    "; of these usable for GEE AR(1) + FG by the paper's filter:", X[, sum(constant & usable)], "\n")
cat("H0, 2,376 design points: mean usable share", round(H0[, mean(pu)], 3), "(Table 4:", round(T4$pu[T4$strategy == "gee_ar1_fg"], 3),
    "); without the constant-outcome datasets", round(H0[, mean(pu_noconst)], 3), "\n")
cat("H0, 2-3 occasions and <= 25 clusters (", nrow(small), "design points): mean usable share", round(small[, mean(pu)], 3),
    "(Table 4:", round(T4$pu_small[T4$strategy == "gee_ar1_fg"], 3), "); without them", round(small[, mean(pu_noconst)], 3), "\n")
cat("largest share of constant-outcome datasets in a design point (H0):", round(H0[, max(const)], 2), "% at",
    H0[which.max(const), paste(dgp, scen, "N", N, "T", ntime)], "\n")
cat("study rows read unchanged (md5 of three files):", identical(md5_before, tools::md5sum(files[c(1, 23760, 47520)])), "\n")
sink()
