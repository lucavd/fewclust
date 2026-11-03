# c1_identity_compare.R -- fewclust 0.4.1, gates A and B on the 60 datasets: compare the files of c1_identity_run.R.
# Usage: Rscript --vanilla c1_identity_compare.R <ref040.rds> <dev041.rds> <out_dir> <constants_040> <constants_041>
args <- commandArgs(TRUE); a <- readRDS(args[1]); b <- readRDS(args[2]); OUT <- args[3]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common_041.R"))
stopifnot(a$version == "0.4.0", b$version == "0.4.1", length(a$res) == length(b$res))
cs <- a$cases; gee_ar1_fg <- which(cs$method == "gee" & cs$corstr == "ar1" & cs$correction == "fg")
rows <- list()
for (s in seq_along(a$res)) for (k in seq_len(nrow(cs))) {
  x <- a$res[[s]]$fits[[k]]; y <- b$res[[s]]$fits[[k]]; g <- b$res[[s]]$fits[[gee_ar1_fg]]
  pred <- predicted_change(x)
  same <- identical(x$refused, y$refused) && (x$refused || (identical(x$method, y$method) && same_table(x$table, y$table)))
  ok_changed <- if (!pred) NA else if (cs$method[k] == "mmm") isTRUE(y$refused) &&
    grepl("did not converge \\(nlm code|no finite estimate", gsub("\\s+", " ", y$message)) else
    if (isTRUE(g$refused)) isTRUE(y$refused) && grepl("did not converge \\(nlm code|no finite estimate", gsub("\\s+", " ", y$message)) else
      !isTRUE(y$refused) && identical(y$method, "gee") && same_table(y$table, g$table)
  ck <- y$convergence
  rp <- rule_path(ck)
  actual <- if (is.null(ck)) NA_character_ else ck$path
  outcome_ok <- is.na(rp) || (rp == actual && (if (rp == "no result") (y$refused && cs$method[k] == "mmm") ||
    (!y$refused && y$method == "gee") || (y$refused && isTRUE(g$refused)) else
    !y$refused || !grepl("convergence check", y$message)))
  fits_same <- if (!is.null(x$convergence) && !is.null(ck)) same_fits(x$convergence, ck) else NA
  rows[[length(rows) + 1L]] <- data.frame(dataset = s, method = cs$method[k], corstr = cs$corstr[k],
    correction = cs$correction[k], coherence = cs$coherence[k], refused_040 = x$refused, refused_041 = y$refused,
    method_040 = if (x$refused) NA else x$method, method_041 = if (y$refused) NA else y$method,
    path_040 = if (is.null(x$convergence)) NA else x$convergence$path, path_041 = actual, rule_path_041 = rp,
    refit_code = if (is.null(ck$constant)) NA else ck$constant$code, predicted_change = pred, identical_040 = same,
    changed_as_expected = ok_changed, rule_ok = outcome_ok, recorded_fits_identical = fits_same)
}
R <- do.call(rbind, rows)
write.csv(R, file.path(OUT, "c1_identity.csv"), row.names = FALSE)
flip_ok <- vapply(seq_along(a$res), function(s) identical(a$res[[s]]$flip, b$res[[s]]$flip), logical(1))
vc_ok <- vapply(seq_along(a$res), function(s) identical(a$res[[s]]$vcov, b$res[[s]]$vcov), logical(1))
const_ok <- unname(tools::md5sum(args[4]) == tools::md5sum(args[5]))
gA <- all(R$identical_040[!R$predicted_change]) && all(R$changed_as_expected[R$predicted_change]) && all(flip_ok) &&
  all(vc_ok) && const_ok && all(R$recorded_fits_identical, na.rm = TRUE)
gB <- all(R$rule_ok)
sink(file.path(OUT, "c1_summary.txt"), split = TRUE)
cat("calls compared:", nrow(R), "(", length(a$res), "datasets x", nrow(cs), "calls); predicted to change:", sum(R$predicted_change),
    "; identical to 0.4.0 among the others:", sum(R$identical_040[!R$predicted_change]), "of", sum(!R$predicted_change), "\n")
cat("  changed as expected (mmm -> fewclust_error naming the reason; auto -> GEE AR(1) FG identical to method = 'gee'):",
    sum(R$changed_as_expected, na.rm = TRUE), "of", sum(R$predicted_change), "\n")
if (any(R$predicted_change)) print(R[R$predicted_change, c("dataset", "method", "coherence", "path_040", "refit_code", "path_041", "method_041")])
bad <- R[(!R$predicted_change & !R$identical_040) | (R$predicted_change & !R$changed_as_expected %in% TRUE) | !R$rule_ok, ]
if (nrow(bad)) { cat("  discrepancies:\n"); print(bad) }
cat("  recorded fits (trend, constant) identical where both versions recorded them:", sum(R$recorded_fits_identical, na.rm = TRUE),
    "of", sum(!is.na(R$recorded_fits_identical)), "\n")
cat("  flip_test identical:", sum(flip_ok), "of", length(flip_ok), "; gee_vcov identical:", sum(vc_ok), "of", length(vc_ok),
    "; R/constants.R md5 identical:", const_ok, "\n")
cat("  rule of D-90 on the package's own recorded fits: consistent in", sum(R$rule_ok), "of", nrow(R), "calls\n")
cat("GATE A (60 datasets):", if (gA) "PASS" else "FAIL", "\n")
cat("GATE B (60 datasets):", if (gB) "PASS" else "FAIL", "\n")
sink()
