# c3_invariance_compare.R -- fewclust 0.3.1, check 7 (part 2): the two files of c3_invariance_run.R must be
# identical() element by element (same bits), apart from the version label.
# Usage: Rscript --vanilla dev/checks_0.3.1/c3_invariance_compare.R <rds_0.3.0> <rds_0.3.1> <out_dir>
args <- commandArgs(TRUE); a <- readRDS(args[1]); b <- readRDS(args[2]); OUT <- args[3]
stopifnot(a$version == "0.3.0", b$version == "0.3.1")
sv <- vapply(names(a$V), function(k) identical(a$V[[k]], b$V[[k]]), logical(1))
sm <- vapply(names(a$M), function(k) identical(a$M[[k]], b$M[[k]]), logical(1))
n_err <- sum(vapply(a$V, function(x) is.character(x) || any(vapply(x, is.character, logical(1))), logical(1)))
sink(file.path(OUT, "c3_summary.txt"), split = TRUE)
cat("gee_vcov (4 estimators x 3 working correlations x scale.fix TRUE/FALSE x 100 datasets):",
    sum(sv), "of", length(sv), "sets identical; sets containing an error message (identical too):", n_err, "\n")
cat("marginalized-model engine with explicit control (40 datasets):", sum(sm), "of", length(sm), "identical\n")
if (!all(sv)) print(head(names(sv)[!sv]))
if (!all(sm)) print(head(names(sm)[!sm]))
cat("versions compared:", a$version, "vs", b$version, "\n")
cat("GATE check 7:", if (all(sv) && all(sm) && identical(names(a$V), names(b$V))) "PASS" else "FAIL", "\n")
sink()
