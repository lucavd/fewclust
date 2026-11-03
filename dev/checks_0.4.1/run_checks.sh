#!/usr/bin/env bash
# run_checks.sh -- fewclust 0.4.1 (dev/aggiornamento_0.4.1.md): install 0.4.0 (reference, from e7bb055, scratch library)
# and 0.4.1 (~/R/fewclust_041_lib), tests, lint, gates A-D, coherence rerun of S9, R CMD check.
# Usage (from ~/MMM-package/fewclust): bash dev/checks_0.4.1/run_checks.sh
# At most 6 concurrent R processes (each mclapply step: 1 parent + 5 workers; steps run one after the other).
set -u
cd "$(dirname "$0")/../.."
S=/home/user/.claude-science-scratch/fc041; OUT=dev/checks_0.4.1/out; mkdir -p "$S/tmp" "$S/lib040" "$OUT"
LIB041=$HOME/R/fewclust_041_lib; LIB040=$S/lib040; mkdir -p "$LIB041"
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TOOLS=/home/user/MMM-package/.Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent TMPDIR="$S/tmp"
rm -rf "$S/src040" && mkdir -p "$S/src040" && git -C /home/user/MMM-package archive e7bb055 fewclust | tar -x -C "$S/src040"
R_LIBS="$PROJ" R CMD INSTALL --no-test-load --library="$LIB040" "$S/src040/fewclust" > "$OUT/install_040.log" 2>&1 || { tail -n 30 "$OUT/install_040.log"; exit 1; }
R_LIBS="$PROJ" R CMD INSTALL --no-test-load --library="$LIB041" . > "$OUT/install.log" 2>&1 || { tail -n 30 "$OUT/install.log"; exit 1; }
echo "installed: $(Rscript --vanilla -e ".libPaths(c('$LIB040','$PROJ')); cat(as.character(packageVersion('fewclust')))") (reference), $(Rscript --vanilla -e ".libPaths(c('$LIB041','$PROJ')); cat(as.character(packageVersion('fewclust')))")"
# gate D: tests and lint
Rscript --vanilla -e ".libPaths(c('$LIB041','$PROJ','$TOOLS')); r <- as.data.frame(testthat::test_local('.', reporter = 'summary', stop_on_failure = FALSE)); cat('\ntests: files', length(unique(r\$file)), '; expectations', sum(r\$nb), '; failed', sum(r\$failed), '; errors', sum(r\$error), '; skipped', sum(r\$skipped), '; warnings', sum(r\$warning), '\n'); write.csv(r[, c('file','test','nb','failed','skipped','error','warning')], '$OUT/tests.csv', row.names = FALSE); cat('GATE D (tests):', if (sum(r\$failed) + sum(r\$error) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/tests.log" 2>&1
tail -n 2 "$OUT/tests.log"
Rscript --vanilla -e ".libPaths(c('$LIB041','$TOOLS','$PROJ')); l <- as.data.frame(lintr::lint_package('.', exclusions = list('R/constants.R', 'tests', 'vignettes', 'dev', 'data-raw'))); write.csv(l, '$OUT/lint.csv', row.names = FALSE); cat('lints in R/:', nrow(l), '\n'); if (nrow(l)) print(table(l\$linter)); cat('GATE D (lint):', if (nrow(l) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/lint.log" 2>&1
tail -n 1 "$OUT/lint.log"
# gates A, B: 60 package datasets
Rscript --vanilla dev/checks_0.4.1/c1_identity_run.R "$LIB040" "$S/datasets.rds" "$S/id_040.rds" > "$OUT/c1_040.log" 2>&1
Rscript --vanilla dev/checks_0.4.1/c1_identity_run.R "$LIB041" "$S/datasets.rds" "$S/id_041.rds" > "$OUT/c1_041.log" 2>&1
Rscript --vanilla dev/checks_0.4.1/c1_identity_compare.R "$S/id_040.rds" "$S/id_041.rds" "$OUT" "$S/src040/fewclust/R/constants.R" R/constants.R > /dev/null 2>&1
tail -n 2 "$OUT/c1_summary.txt"
# gates A, B, C: study sample (drawn once by c2_study_sample.R; redrawn only if missing)
[ -f "$OUT/c2_sample.csv" ] && [ -f "$S/c2_md5_sources_before.csv" ] || Rscript --vanilla dev/checks_0.4.1/c2_study_sample.R "$OUT" "$S" > "$OUT/c2_sample.log" 2>&1
Rscript --vanilla dev/checks_0.4.1/c2_study_run.R "$LIB040" "$OUT/c2_sample.csv" "$S/st_040.rds" > "$OUT/c2_040.log" 2>&1
Rscript --vanilla dev/checks_0.4.1/c2_study_run.R "$LIB041" "$OUT/c2_sample.csv" "$S/st_041.rds" > "$OUT/c2_041.log" 2>&1
Rscript --vanilla dev/checks_0.4.1/c2_study_compare.R "$S/st_040.rds" "$S/st_041.rds" "$OUT/c2_sample.csv" "$S/c2_md5_sources_before.csv" "$OUT" > /dev/null 2>&1
tail -n 3 "$OUT/c2_summary.txt"
# prediction 4: the S9 coherence comparison with 0.4.1
Rscript --vanilla dev/checks_0.4.1/c3_coherence.R "$LIB041" "$OUT" > "$OUT/c3.log" 2>&1
grep -E "^compared|^description" "$OUT/c3_summary.txt"
# gate D: R CMD check
bash dev/checks_0.4.1/c7_rcmdcheck.sh "$LIB041" "$OUT" "$S/tmp"
