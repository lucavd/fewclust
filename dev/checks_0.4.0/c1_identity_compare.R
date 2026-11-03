# c1_identity_compare.R -- fewclust 0.4.0, checks 1 and 2: compare the two files of c1_identity_run.R.
# Usage: Rscript --vanilla c1_identity_compare.R <ref.rds> <dev.rds> <out_dir> <constants_ref> <constants_dev>
args <- commandArgs(TRUE); a <- readRDS(args[1]); b <- readRDS(args[2]); OUT <- args[3]
stopifnot(a$version == "0.3.1", b$version == "0.4.0", length(a$res) == length(b$res))
cols <- c("term", "estimate", "std.error", "statistic", "df", "p.value")
rows <- list(); n_fit <- 0; n_ref <- 0; bad <- character(0)
for (s in seq_along(a$res)) for (k in seq_len(nrow(a$cases))) {
  x <- a$res[[s]]$fits[[k]]; y <- b$res[[s]]$fits[[k]]
  ok <- if (x$refused) { n_ref <- n_ref + 1; isTRUE(y$refused) } else {
    n_fit <- n_fit + 1
    !isTRUE(y$refused) && identical(x$method, y$method) &&
      identical(unname(as.list(x$table[, cols])), unname(as.list(y$table[, cols])))
  }
  if (!ok) bad <- c(bad, sprintf("dataset %d, %s/%s/%s", s, a$cases$method[k], a$cases$corstr[k], a$cases$correction[k]))
}
flip_ok <- vapply(seq_along(a$res), function(s) identical(a$res[[s]]$flip, b$res[[s]]$flip), logical(1))
vc_ok <- vapply(seq_along(a$res), function(s) identical(a$res[[s]]$vcov, b$res[[s]]$vcov), logical(1))
const_ok <- unname(tools::md5sum(args[4]) == tools::md5sum(args[5]))
sink(file.path(OUT, "c1_summary.txt"), split = TRUE)
cat("check 1: fits compared with 0.3.1:", n_fit, "; refusals in 0.3.1 (0.4.0 must error on each):", n_ref,
    "; mismatches:", length(bad), "\n")
if (length(bad)) print(head(bad, 20))
cat("  flip_test p-value and reason identical:", sum(flip_ok), "of", length(flip_ok), "\n")
cat("check 2: gee_vcov identical:", sum(vc_ok), "of", length(vc_ok), "datasets (3 corstr x 2 scale.fix x 4 types each);",
    "R/constants.R md5 identical:", const_ok, "\n")
cat("GATE check 1:", if (!length(bad) && all(flip_ok)) "PASS" else "FAIL", "\n")
cat("GATE check 2:", if (all(vc_ok) && const_ok) "PASS" else "FAIL", "\n")
sink()
