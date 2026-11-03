#!/usr/bin/env bash
# c7_rcmdcheck.sh -- fewclust 0.3.0, check 4: R CMD build + R CMD check --as-cran, no errors and no warnings.
# Usage (from the package root): bash dev/checks_0.3.0/c7_rcmdcheck.sh <dev_lib> <out_dir>
set -u
DEV=$(cd "$1" && pwd); OUT=$(mkdir -p "$2" && cd "$2" && pwd)
PROJ=/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
TOOLS=/home/user/MMM-package/.Rlib
export R_LIBS="$DEV:$PROJ:$TOOLS" R_LIBS_USER=/nonexistent R_LIBS_SITE=/nonexistent
export _R_CHECK_CRAN_INCOMING_REMOTE_=false OMP_NUM_THREADS=1 MAKEFLAGS=-j4
# ~/bin holds a `tidy` wrapper that points to a missing temporary file of another session (it aborts the check of the
# HTML manual); ~/bin also holds pdflatex, so it stays in PATH and only R's tidy command is pointed elsewhere
export R_TIDYCMD=tidy-not-installed
PKG=$(pwd); TMP=$(mktemp -d); cd "$TMP"
R CMD build "$PKG" > "$OUT/c7_build.log" 2>&1
TGZ=$(ls fewclust_*.tar.gz | head -n 1)
R CMD check --as-cran "$TGZ" > "$OUT/c7_check.log" 2>&1
tail -n 40 "$OUT/c7_check.log" > "$OUT/c7_summary.txt"
grep -E "^Status:|WARNING|ERROR|NOTE" "$OUT/c7_check.log" >> "$OUT/c7_summary.txt" || true
if grep -q "^Status: OK$\|^Status: [0-9]* NOTEs*$" "$OUT/c7_check.log"; then echo "GATE check 4: PASS" >> "$OUT/c7_summary.txt"; else echo "GATE check 4: FAIL" >> "$OUT/c7_summary.txt"; fi
cp -r fewclust.Rcheck/00check.log "$OUT/c7_00check.log" 2>/dev/null || true
cat "$OUT/c7_summary.txt"
rm -rf "$TMP"
