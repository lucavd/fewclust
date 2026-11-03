#!/usr/bin/env bash
# run_checks.sh -- fewclust 0.4.0: install, tests, identity with 0.3.1, lint, R CMD check.
# Usage (from ~/MMM-package/fewclust): bash dev/checks_0.4.0/run_checks.sh
# Reference 0.3.1 must be installed in ~/MMM-package/.tmp040/lib031 (from commit 626393b). At most 4 cores.
set -u
cd "$(dirname "$0")/../.."
DEV=$HOME/R/fewclust_dev_lib; OUT=dev/checks_0.4.0/out; mkdir -p "$DEV" "$OUT"
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TOOLS=/home/user/MMM-package/.Rlib; T=/home/user/MMM-package/.tmp040; mkdir -p "$T/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent TMPDIR="$T/tmp"
R_LIBS="$DEV:$PROJ:$TOOLS" R CMD INSTALL --no-test-load --library="$DEV" . > "$OUT/install.log" 2>&1 || { tail -n 30 "$OUT/install.log"; exit 1; }
echo "installed: $(Rscript --vanilla -e ".libPaths(c('$DEV','$PROJ')); cat(as.character(packageVersion('fewclust')))")"
Rscript --vanilla -e ".libPaths(c('$DEV','$PROJ','$TOOLS')); r <- as.data.frame(testthat::test_local('.', reporter = 'summary', stop_on_failure = FALSE)); cat('\ntests: files', length(unique(r\$file)), '; expectations', sum(r\$nb), '; failed', sum(r\$failed), '; errors', sum(r\$error), '; skipped', sum(r\$skipped), '; warnings', sum(r\$warning), '\n'); write.csv(r[, c('file','test','nb','failed','skipped','error','warning')], '$OUT/tests.csv', row.names = FALSE); cat('GATE check 3-5 (tests):', if (sum(r\$failed) + sum(r\$error) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/tests.log" 2>&1
tail -n 2 "$OUT/tests.log"
Rscript --vanilla dev/checks_0.4.0/c1_identity_run.R "$T/lib031" "$T/datasets.rds" "$T/id_031.rds" > "$OUT/c1_031.log" 2>&1
Rscript --vanilla dev/checks_0.4.0/c1_identity_run.R "$DEV" "$T/datasets.rds" "$T/id_040.rds" > "$OUT/c1_040.log" 2>&1
Rscript --vanilla dev/checks_0.4.0/c1_identity_compare.R "$T/id_031.rds" "$T/id_040.rds" "$OUT" "$T/src031/fewclust/R/constants.R" R/constants.R > /dev/null 2>&1
tail -n 2 "$OUT/c1_summary.txt"
Rscript --vanilla -e ".libPaths(c('$DEV','$TOOLS','$PROJ')); l <- as.data.frame(lintr::lint_package('.', exclusions = list('R/constants.R', 'tests', 'vignettes', 'dev', 'data-raw'))); write.csv(l, '$OUT/lint.csv', row.names = FALSE); cat('lints in R/:', nrow(l), '\n'); if (nrow(l)) print(table(l\$linter)); cat('GATE check 7:', if (nrow(l) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/lint.log" 2>&1
tail -n 1 "$OUT/lint.log"
bash dev/checks_0.4.0/c7_rcmdcheck.sh "$DEV" "$OUT" "$T/tmp"
