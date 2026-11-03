#!/usr/bin/env bash
# run_checks.sh -- install the working copy into the development library and run every check of the 0.3.0 update.
# Usage (from ~/MMM-package/fewclust): bash dev/checks_0.3.0/run_checks.sh [which]   which = all | tests | study | fast
# The development library is ~/R/fewclust_dev_lib (never ~/R/s7d_lib, read by the paper pipeline). Uses at most 4 cores.
set -u
cd "$(dirname "$0")/../.."
WHICH=${1:-all}
DEV=$HOME/R/fewclust_dev_lib; OUT=dev/checks_0.3.0/out; mkdir -p "$DEV" "$OUT"
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TOOLS=/home/user/MMM-package/.Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent
MMM_HEAD=$(git -C ~/MMM rev-parse HEAD); MMM_STATUS=$(git -C ~/MMM status --porcelain | md5sum)
R_LIBS="$DEV:$PROJ:$TOOLS" R CMD INSTALL --no-test-load --library="$DEV" . > "$OUT/install.log" 2>&1 || { tail -n 30 "$OUT/install.log"; exit 1; }
echo "installed: $(Rscript --vanilla -e ".libPaths(c('$DEV','$PROJ')); cat(as.character(packageVersion('fewclust')))")"
if [ "$WHICH" = all ] || [ "$WHICH" = tests ] || [ "$WHICH" = fast ]; then
  Rscript --vanilla -e ".libPaths(c('$DEV','$PROJ','$TOOLS')); options(testthat.progress.max_fails = 1000); r <- as.data.frame(testthat::test_local('.', reporter = 'summary', stop_on_failure = FALSE)); cat('\ntests: files', length(unique(r\$file)), '; expectations', sum(r\$nb), '; failed', sum(r\$failed), '; errors', sum(r\$error), '; skipped', sum(r\$skipped), '\n'); write.csv(r[, c('file','test','nb','failed','skipped','error')], '$OUT/tests.csv', row.names = FALSE)" > "$OUT/tests.log" 2>&1
  tail -n 3 "$OUT/tests.log"
fi
if [ "$WHICH" = all ] || [ "$WHICH" = fast ]; then
  Rscript --vanilla dev/checks_0.3.0/c2_rule_fig1.R "$DEV" "$OUT" > /dev/null 2>&1; tail -n 1 "$OUT/c2_summary.txt"
  Rscript --vanilla dev/checks_0.3.0/c4_paper_condition.R "$DEV" "$OUT" > /dev/null 2>&1; tail -n 1 "$OUT/c4_summary.txt"
  Rscript --vanilla dev/checks_0.3.0/c5_constants.R "$DEV" "$OUT" > /dev/null 2>&1; tail -n 1 "$OUT/c5_summary.txt"
  Rscript --vanilla dev/checks_0.3.0/c3_invariance_run.R "$HOME/R/s7d_lib" /tmp/fc_inv_020.rds > /dev/null 2>&1
  Rscript --vanilla dev/checks_0.3.0/c3_invariance_run.R "$DEV" /tmp/fc_inv_030.rds > /dev/null 2>&1
  Rscript --vanilla dev/checks_0.3.0/c3_invariance_compare.R /tmp/fc_inv_020.rds /tmp/fc_inv_030.rds "$OUT" > /dev/null 2>&1; tail -n 1 "$OUT/c3_summary.txt"
fi
if [ "$WHICH" = all ] || [ "$WHICH" = study ]; then
  Rscript --vanilla dev/checks_0.3.0/c1_study_repro.R "$DEV" "$OUT" > "$OUT/c1.log" 2>&1; tail -n 2 "$OUT/c1_summary.txt"
  Rscript --vanilla dev/checks_0.3.0/c6_coherence.R "$DEV" "$OUT" > "$OUT/c6.log" 2>&1; grep "P1\|relevance" "$OUT/c6_summary.txt"
fi
if [ "$WHICH" = all ]; then bash dev/checks_0.3.0/c7_rcmdcheck.sh "$DEV" "$OUT" | tail -n 1; fi
echo "~/MMM untouched: $( [ "$(git -C ~/MMM rev-parse HEAD)" = "$MMM_HEAD" ] && [ "$(git -C ~/MMM status --porcelain | md5sum)" = "$MMM_STATUS" ] && echo yes || echo NO )"
