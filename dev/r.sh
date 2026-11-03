#!/usr/bin/env bash
# Run R code with the package development libraries, bypassing the repo .Rprofile.
# Usage: dev/r.sh 'any R code'
LIBS='.libPaths(c("/home/user/MMM-package/.Rlib","/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))'
OPTS='options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"), HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])), Ncpus = 8)'
exec Rscript --vanilla -e "$LIBS" -e "$OPTS" -e "$1"
