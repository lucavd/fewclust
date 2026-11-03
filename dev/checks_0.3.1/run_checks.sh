#!/usr/bin/env bash
# run_checks.sh -- install the working copy into the development library and run every check of the 0.3.1 update.
# Usage (from ~/MMM-package/fewclust): bash dev/checks_0.3.1/run_checks.sh [which]   which = all | tests | study | fast | ref
#   ref: install fewclust 0.3.0 (commit 22e618e) into ~/MMM-package/.tmp031/lib030, the reference of check 7.
# The development library is ~/R/fewclust_dev_lib (never ~/R/s7d_lib, read by the paper pipeline). Uses at most 4 cores.
# Temporary files live in ~/MMM-package/.tmp031 (never /tmp): the update writes only under ~/MMM-package.
set -u
cd "$(dirname "$0")/../.."
WHICH=${1:-all}
DEV=$HOME/R/fewclust_dev_lib; OUT=dev/checks_0.3.1/out; mkdir -p "$DEV" "$OUT"
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TOOLS=/home/user/MMM-package/.Rlib
TROOT=/home/user/MMM-package/.tmp031; mkdir -p "$TROOT/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent TMPDIR="$TROOT/tmp"
INPUTS="/home/user/MMM/tables/V5_S7f_table4.csv /home/user/MMM/tables/V5_S7f_h2h_boot.csv /home/user/MMM/tables/V5_S7f_rule_dp.csv /home/user/MMM/tables/V5_S7f_cp_S7d_dp_rates.csv"
MD5_IN=$(md5sum $INPUTS)
echo "inputs at start:"; echo "$MD5_IN"
echo "~/MMM at start (informative; another session works there): $(git -C ~/MMM rev-parse --short HEAD), $(git -C ~/MMM status --porcelain | wc -l) changed paths"
if [ "$WHICH" = ref ]; then
  rm -rf "$TROOT/src030" "$TROOT/lib030"; mkdir -p "$TROOT/src030" "$TROOT/lib030"
  git -C /home/user/MMM-package archive 22e618e fewclust | tar -x -C "$TROOT/src030"
  R_LIBS="$PROJ:$TOOLS" R CMD INSTALL --no-test-load --library="$TROOT/lib030" "$TROOT/src030/fewclust" > "$OUT/install_030.log" 2>&1 || { tail -n 30 "$OUT/install_030.log"; exit 1; }
  echo "reference installed: $(Rscript --vanilla -e ".libPaths(c('$TROOT/lib030','$PROJ')); cat(as.character(packageVersion('fewclust')))")"
  exit 0
fi
R_LIBS="$DEV:$PROJ:$TOOLS" R CMD INSTALL --no-test-load --library="$DEV" . > "$OUT/install.log" 2>&1 || { tail -n 30 "$OUT/install.log"; exit 1; }
echo "installed: $(Rscript --vanilla -e ".libPaths(c('$DEV','$PROJ')); cat(as.character(packageVersion('fewclust')))")"
if [ "$WHICH" = all ] || [ "$WHICH" = tests ] || [ "$WHICH" = fast ]; then
  Rscript --vanilla -e ".libPaths(c('$DEV','$PROJ','$TOOLS')); options(testthat.progress.max_fails = 1000); r <- as.data.frame(testthat::test_local('.', reporter = 'summary', stop_on_failure = FALSE)); cat('\ntests: files', length(unique(r\$file)), '; expectations', sum(r\$nb), '; failed', sum(r\$failed), '; errors', sum(r\$error), '; skipped', sum(r\$skipped), '\n'); write.csv(r[, c('file','test','nb','failed','skipped','error')], '$OUT/tests.csv', row.names = FALSE)" > "$OUT/tests.log" 2>&1
  tail -n 3 "$OUT/tests.log"
fi
if [ "$WHICH" = all ] || [ "$WHICH" = fast ]; then
  Rscript --vanilla dev/checks_0.3.1/c2_rule_fig1.R "$DEV" "$OUT" > "$OUT/c2.log" 2>&1; tail -n 1 "$OUT/c2_summary.txt"
  Rscript --vanilla dev/checks_0.3.1/c4_paper_condition.R "$DEV" "$OUT" > "$OUT/c4.log" 2>&1; tail -n 1 "$OUT/c4_summary.txt"
  Rscript --vanilla dev/checks_0.3.1/c5_constants.R "$DEV" "$OUT" > "$OUT/c5.log" 2>&1; tail -n 1 "$OUT/c5_summary.txt"
  Rscript --vanilla dev/checks_0.3.1/c8_constant_outcome.R "$DEV" "$OUT" > "$OUT/c8.log" 2>&1; tail -n 1 "$OUT/c8_summary.txt"
  Rscript --vanilla dev/checks_0.3.1/c3_invariance_run.R "$TROOT/lib030" "$TROOT/fc_inv_030.rds" > "$OUT/c3_030.log" 2>&1
  Rscript --vanilla dev/checks_0.3.1/c3_invariance_run.R "$DEV" "$TROOT/fc_inv_031.rds" > "$OUT/c3_031.log" 2>&1
  Rscript --vanilla dev/checks_0.3.1/c3_invariance_compare.R "$TROOT/fc_inv_030.rds" "$TROOT/fc_inv_031.rds" "$OUT" > /dev/null 2>&1; tail -n 1 "$OUT/c3_summary.txt"
fi
if [ "$WHICH" = all ] || [ "$WHICH" = study ]; then
  Rscript --vanilla dev/checks_0.3.1/c1_study_repro.R "$DEV" "$OUT" > "$OUT/c1.log" 2>&1; tail -n 2 "$OUT/c1_summary.txt"
  Rscript --vanilla dev/checks_0.3.1/c9_mmm_repro.R "$DEV" "$OUT" > "$OUT/c9.log" 2>&1; tail -n 1 "$OUT/c9_summary.txt"
fi
if [ "$WHICH" = all ]; then bash dev/checks_0.3.1/c7_rcmdcheck.sh "$DEV" "$OUT" | tail -n 1; fi
echo "inputs unchanged during the checks: $( [ "$(md5sum $INPUTS)" = "$MD5_IN" ] && echo yes || echo NO )"
echo "~/MMM at end (informative): $(git -C ~/MMM rev-parse --short HEAD), $(git -C ~/MMM status --porcelain | wc -l) changed paths"
