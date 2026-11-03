#!/usr/bin/env bash
# c7_rcmdcheck.sh -- fewclust 0.4.0, check 6: R CMD build + R CMD check --as-cran, no errors and no warnings.
set -u
DEV=$(cd "$1" && pwd); OUT=$(mkdir -p "$2" && cd "$2" && pwd); TMPROOT=$3
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TOOLS=/home/user/MMM-package/.Rlib
export R_LIBS="$DEV:$PROJ:$TOOLS" R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent
export _R_CHECK_CRAN_INCOMING_REMOTE_=false OMP_NUM_THREADS=1 MAKEFLAGS=-j4
# ~/bin/tidy is a broken wrapper of another session; ~/bin also holds pdflatex, so only R's tidy is redirected
export R_TIDYCMD=tidy-not-installed
PKG=$(pwd); TMP=$(mktemp -d -p "$TMPROOT"); cd "$TMP"
R CMD build "$PKG" > "$OUT/c7_build.log" 2>&1
TGZ=$(ls fewclust_*.tar.gz | head -n 1)
R CMD check --as-cran "$TGZ" > "$OUT/c7_check.log" 2>&1
grep -E "^Status:|WARNING|ERROR|NOTE" "$OUT/c7_check.log" > "$OUT/c7_summary.txt" || true
if grep -q "^Status: OK$\|^Status: [0-9]* NOTEs*$" "$OUT/c7_check.log"; then echo "GATE check 6: PASS" >> "$OUT/c7_summary.txt"; else echo "GATE check 6: FAIL" >> "$OUT/c7_summary.txt"; fi
tail -n 1 "$OUT/c7_summary.txt"
rm -rf "$TMP"
