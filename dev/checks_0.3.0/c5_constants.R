# c5_constants.R -- fewclust 0.3.0, check 8: every count and share the package quotes equals Table 4 of the paper
# (branch morris-v4, tables/V5_S7f_table4.csv, read here, not retyped), and no text of the package still quotes the
# 1,896 design points of manuscript v4.
# Usage (from the package root): Rscript --vanilla dev/checks_0.3.0/c5_constants.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(fewclust)); stopifnot(as.character(packageVersion("fewclust")) == "0.3.0")
T4 <- read.csv(path.expand("~/MMM/tables/V5_S7f_table4.csv"))
st <- fewclust:::.fc_study
P <- st$table4
cols <- c("hi_all", "lo_all", "lo_N10", "n_N10", "pu", "pu_small", "n_dp")
both <- intersect(T4$strategy, rownames(P))
cmp <- do.call(rbind, lapply(both, function(s) data.frame(strategy = s, field = cols,
  table4 = unlist(T4[T4$strategy == s, cols]), package = unlist(P[s, cols]))))
cmp$equal <- abs(cmp$table4 - cmp$package) <= 1e-9 * pmax(1, abs(cmp$table4))
other <- c(design_points = st$design_points == unique(T4$n_dp), all_strategies = setequal(both, T4$strategy))
txt <- c(list.files("R", full.names = TRUE), list.files("man", full.names = TRUE), list.files("vignettes", "Rmd$", full.names = TRUE),
         "README.md", "DESCRIPTION", "inst/CITATION")
stale <- do.call(rbind, lapply(txt, function(f) { l <- readLines(f, warn = FALSE); k <- grep("1,896|1896|manuscript v4", l)
  if (length(k)) data.frame(file = f, line = k, text = substr(l[k], 1, 120)) else NULL }))
write.csv(cmp, file.path(OUT, "c5_constants.csv"), row.names = FALSE)
sink(file.path(OUT, "c5_summary.txt"), split = TRUE)
cat("Table 4 values compared:", nrow(cmp), "; equal:", sum(cmp$equal), "; strategies:", length(both), "of", nrow(T4), "\n")
print(other); if (!all(cmp$equal)) print(cmp[!cmp$equal, ])
cat("lines quoting the v4 study (1,896 design points):", if (is.null(stale)) 0L else nrow(stale), "\n"); if (!is.null(stale)) print(stale)
cat("version in DESCRIPTION:", read.dcf("DESCRIPTION")[1, "Version"], "; NEWS opens with 0.3.0:", grepl("0.3.0", readLines("NEWS.md")[1]), "\n")
ok <- all(cmp$equal) && all(other) && is.null(stale) && read.dcf("DESCRIPTION")[1, "Version"] == "0.3.0" && grepl("0.3.0", readLines("NEWS.md")[1])
cat("GATE check 8:", if (ok) "PASS" else "FAIL", "\n"); sink()
