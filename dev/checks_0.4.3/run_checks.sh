#!/usr/bin/env bash
# run_checks.sh -- fewclust 0.4.3 (dev/aggiornamento_0.4.3.md): regenerate R/constants.R, roxygen, gates 1-5.
# Usage (from ~/MMM-package/fewclust): bash dev/checks_0.4.3/run_checks.sh
set -u
cd "$(dirname "$0")/../.."
OUT=dev/checks_0.4.3/out; S=/home/user/.claude-science-scratch/fc043; rm -rf "$S"; mkdir -p "$S/tmp" "$S/lib" "$OUT"
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu; TOOLS=/home/user/MMM-package/.Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent TMPDIR="$S/tmp" PATH="$HOME/bin:$PATH"
git show f05c914:fewclust/R/constants.R > "$S/constants_042.R"
Rscript --vanilla dev/make_constants.R > "$OUT/make_constants.log" 2>&1; tail -n 1 "$OUT/make_constants.log"
diff "$S/constants_042.R" R/constants.R > "$OUT/constants_vs_042.diff"
NV=$(grep -c '^[<>]' "$OUT/constants_vs_042.diff"); NH=$(grep '^[<>]' "$OUT/constants_vs_042.diff" | grep -cE '^[<>] # lucavd/MMM \(|^[<>] #   [0-9a-f]{32}  V5_S7f_numbers.csv$')
echo "constants vs 0.4.2: changed lines $NV, of which header lines (commit, registry md5) $NH" | tee "$OUT/c1_summary.txt"
Rscript --vanilla -e ".libPaths(c('$PROJ','$TOOLS')); roxygen2::roxygenise('.')" > "$OUT/roxygen.log" 2>&1
R_LIBS="$PROJ" R CMD INSTALL --no-test-load --library="$S/lib" . > "$OUT/install.log" 2>&1 || { tail -n 30 "$OUT/install.log"; exit 1; }
mkdir -p "$S/c1_042"; Rscript --vanilla dev/checks_0.4.2/c1_constants.R "$S/lib" "$S/constants_042.R" "$S/c1_042" > "$OUT/c1_constants.log" 2>&1
cat "$S/c1_042/c1_summary.txt" >> "$OUT/c1_summary.txt" 2>/dev/null
OK1=$([ "$NV" = "$NH" ] && grep -q "GATE 1: PASS" "$OUT/c1_summary.txt" && echo PASS || echo FAIL); echo "GATE 1 (constants): $OK1" | tee -a "$OUT/c1_summary.txt"
CH=$(git diff --name-only f05c914 -- R | sort | tr '\n' ' '); echo "R/ files changed since 0.4.2: $CH" | tee "$OUT/c2_summary.txt"
if [ "$CH" = "fewclust/R/constants.R fewclust/R/methods.R " ]; then echo "GATE 2: PASS" | tee -a "$OUT/c2_summary.txt"; else echo "GATE 2: FAIL" | tee -a "$OUT/c2_summary.txt"; fi
Rscript --vanilla -e ".libPaths(c('$S/lib','$PROJ','$TOOLS')); r <- as.data.frame(testthat::test_local('.', reporter = 'summary', stop_on_failure = FALSE)); cat('\ntests: files', length(unique(r\$file)), '; expectations', sum(r\$nb), '; failed', sum(r\$failed), '; errors', sum(r\$error), '; skipped', sum(r\$skipped), '; warnings', sum(r\$warning), '\n'); write.csv(r[, c('file','test','nb','failed','skipped','error','warning')], '$OUT/tests.csv', row.names = FALSE); cat('GATE 3 (tests):', if (sum(r\$failed) + sum(r\$error) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/tests.log" 2>&1
tail -n 2 "$OUT/tests.log"
Rscript --vanilla -e ".libPaths(c('$S/lib','$TOOLS','$PROJ')); l <- as.data.frame(lintr::lint_package('.', exclusions = list('R/constants.R', 'tests', 'vignettes', 'dev', 'data-raw'))); write.csv(l, '$OUT/lint.csv', row.names = FALSE); cat('lints in R/:', nrow(l), '\n'); cat('GATE 3 (lint):', if (nrow(l) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/lint.log" 2>&1
tail -n 1 "$OUT/lint.log"
Rscript --vanilla -e ".libPaths(c('$PROJ','$TOOLS')); s <- spelling::spell_check_package('.'); print(s); cat('words not in WORDLIST:', nrow(s), '\n'); cat('GATE 5 (spelling):', if (nrow(s) == 0) 'PASS' else 'FAIL', '\n')" > "$OUT/spelling.log" 2>&1
tail -n 1 "$OUT/spelling.log"
PKG=$(pwd); cd "$S"
R_LIBS="$PROJ:$TOOLS" R CMD build "$PKG" > "$PKG/$OUT/build.log" 2>&1
TGZ=$(ls fewclust_*.tar.gz | head -n 1); ls -l "$TGZ" | awk '{print "tarball:", $5, "bytes", $9}' > "$PKG/$OUT/tarball.txt"
export R_LIBS="$PROJ:$TOOLS" R_TIDYCMD=tidy-not-installed _R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE _R_CHECK_FORCE_SUGGESTS_=TRUE MAKEFLAGS=-j1
R CMD check --as-cran --timings "$TGZ" > "$PKG/$OUT/check_console.log" 2>&1
cp fewclust.Rcheck/00check.log "$PKG/$OUT/00check.log"; cp fewclust.Rcheck/fewclust-Ex.timings "$PKG/$OUT/examples_timings.txt"
cd "$PKG"
ST=$(grep '^Status:' "$OUT/00check.log"); NE=$(grep -c '\.\.\. ERROR' "$OUT/00check.log"); NW=$(grep -c '\.\.\. WARNING' "$OUT/00check.log")
{ grep -E '(NOTE|WARNING|ERROR)$' "$OUT/00check.log"; echo "$ST"; echo "GATE 4 (R CMD check --as-cran, incoming remote): $([ "$NE" = 0 ] && [ "$NW" = 0 ] && echo PASS || echo FAIL)"; } > "$OUT/c4_summary.txt"
cat "$OUT/c4_summary.txt"
