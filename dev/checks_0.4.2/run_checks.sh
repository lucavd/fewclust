#!/usr/bin/env bash
# run_checks.sh -- fewclust 0.4.2 (dev/aggiornamento_0.4.2.md): regenerate R/constants.R, install, gates 1-4.
# Usage (from ~/MMM-package/fewclust): bash dev/checks_0.4.2/run_checks.sh <old_constants.R (0.4.1)>
set -u
cd "$(dirname "$0")/../.."
OLD=$1; S=/home/user/.claude-science-scratch/fc041; OUT=dev/checks_0.4.2/out; mkdir -p "$S/tmp" "$OUT"
LIB=$HOME/R/fewclust_042_lib; mkdir -p "$LIB"
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu; TOOLS=/home/user/MMM-package/.Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent TMPDIR="$S/tmp"
Rscript --vanilla dev/make_constants.R > "$OUT/make_constants.log" 2>&1; tail -n 1 "$OUT/make_constants.log"
R_LIBS="$PROJ" R CMD INSTALL --no-test-load --library="$LIB" . > "$OUT/install.log" 2>&1 || { tail -n 30 "$OUT/install.log"; exit 1; }
Rscript --vanilla dev/checks_0.4.2/c1_constants.R "$LIB" "$OLD" "$OUT" > "$OUT/c1.log" 2>&1; tail -n 1 "$OUT/c1_summary.txt"
CH=$(git diff --name-only 45c7946 -- R | tr '\n' ' '); echo "R/ files changed since 0.4.1: $CH"
if [ "$CH" = "fewclust/R/constants.R " ]; then echo "GATE 2: PASS"; else echo "GATE 2: FAIL"; fi
Rscript --vanilla -e ".libPaths(c('$LIB','$PROJ','$TOOLS')); r <- as.data.frame(testthat::test_local('.', reporter = 'summary', stop_on_failure = FALSE)); cat('\ntests: files', length(unique(r\$file)), '; expectations', sum(r\$nb), '; failed', sum(r\$failed), '; errors', sum(r\$error), '; skipped', sum(r\$skipped), '; warnings', sum(r\$warning), '\n'); write.csv(r[, c('file','test','nb','failed','skipped','error','warning')], '$OUT/tests.csv', row.names = FALSE); cat('GATE 3 (tests):', if (sum(r\$failed) + sum(r\$error) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/tests.log" 2>&1
tail -n 2 "$OUT/tests.log"
Rscript --vanilla -e ".libPaths(c('$LIB','$TOOLS','$PROJ')); l <- as.data.frame(lintr::lint_package('.', exclusions = list('R/constants.R', 'tests', 'vignettes', 'dev', 'data-raw'))); write.csv(l, '$OUT/lint.csv', row.names = FALSE); cat('lints in R/:', nrow(l), '\n'); cat('GATE 3 (lint):', if (nrow(l) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/lint.log" 2>&1
tail -n 1 "$OUT/lint.log"
bash dev/checks_0.4.1/c7_rcmdcheck.sh "$LIB" "$OUT" "$S/tmp"
# the rendered vignette, for reading
cp "$OUT/../../checks_0.4.2/out/c7_check.log" /dev/null 2>&1 || true
