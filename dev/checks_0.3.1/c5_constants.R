# c5_constants.R -- fewclust 0.3.1, check 9: every value the package quotes equals the files of the corrected study
# (branch morris-v4, working copy of step S9, read here, not retyped): the per-strategy table (tables/V5_S7f_table4.csv,
# Table 5 of manuscript v6), the chosen strategy, the design points where GEE AR(1) + Fay-Graubard is clearly above 7.5%
# (tables/V5_S7f_rule_dp.csv with counts from tables/V5_S7f_cp_S7d_dp_rates.csv, t reference) and the Figure 1a map
# (tables/V5_S7f_h2h_boot.csv); the md5 of the four files are those of dev/aggiornamento_0.3.1.md; no text of the package
# still carries the old recommendation; DESCRIPTION 0.3.1 and NEWS opening with 0.3.1.
# Usage (from the package root): Rscript --vanilla dev/checks_0.3.1/c5_constants.R <dev_lib> <out_dir>
args <- commandArgs(TRUE); DEV <- normalizePath(args[1]); OUT <- args[2]
.libPaths(c(DEV, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(fewclust)); stopifnot(as.character(packageVersion("fewclust")) == "0.3.1")
TB <- path.expand("~/MMM/tables")
MD5 <- c(V5_S7f_table4.csv = "bae9d66a6dda8790e8cdd1537a062888", V5_S7f_h2h_boot.csv = "14eeef1dbc4c46cd9692b2420331d3bf",
         V5_S7f_rule_dp.csv = "494a22502e437a372f7210a7bb3f39cf", V5_S7f_cp_S7d_dp_rates.csv = "e3b6711c6f423493081ce98f422444d9")
md5_now <- tools::md5sum(file.path(TB, names(MD5))); names(md5_now) <- names(MD5)
T4 <- read.csv(file.path(TB, "V5_S7f_table4.csv"))
st <- fewclust:::.fc_study
eq <- function(a, b) if (is.numeric(a)) abs(a - b) <= 1e-9 * pmax(1, abs(a)) else a == b
# 1. the per-strategy table
P <- st$table4
cols <- setdiff(names(T4)[vapply(T4, is.numeric, logical(1))], character(0))
both <- intersect(T4$strategy, rownames(P))
cmp <- do.call(rbind, lapply(both, function(s) data.frame(what = "table4", key = s, field = cols,
  file = unlist(T4[T4$strategy == s, cols]), package = unlist(P[s, cols]))))
# 2. the GEE AR(1) + Fay-Graubard points clearly above 7.5%
D <- read.csv(file.path(TB, "V5_S7f_rule_dp.csv"))
K <- read.csv(file.path(TB, "V5_S7f_cp_S7d_dp_rates.csv"))
hi <- D[D$strategy == "gee_ar1_fg" & D$hi, c("dgp", "scen", "N", "ntime", "pu")]
hi <- merge(hi, K[K$strategy == "gee_ar1_fg" & K$ref == "t", c("dgp", "scen", "N", "ntime", "n", "nd", "kk", "rate")])
hi <- hi[order(hi$dgp), ]; A <- st$ar1_fg_hi; A <- A[order(A$dgp), ]
cmp2 <- do.call(rbind, lapply(c("N", "ntime", "pu", "n", "nd", "kk", "rate"), function(f)
  data.frame(what = "ar1_fg_hi", key = paste(hi$dgp, hi$scen), field = f, file = hi[[f]], package = A[[f]][seq_len(nrow(hi))])))
same_points <- nrow(hi) == nrow(A) && all(hi$dgp == A$dgp) && all(hi$scen == A$scen) && nrow(hi) == T4$hi_all[T4$strategy == "gee_ar1_fg"]
# 3. the Figure 1a map
H <- read.csv(file.path(TB, "V5_S7f_h2h_boot.csv"))
M <- st$h2h
cmp3 <- do.call(rbind, lapply(M$region, function(r) {
  g <- function(q, col = "est") H[H$region == r & H$q == q, col]
  data.frame(what = "h2h", key = r, field = c("n_dp", "mad_gee", "mad_mmm", "diff", "margin"),
             file = c(g("mad_gee", "n_dp"), g("mad_gee"), g("mad_mmm"), g("mad_gee_minus_mmm"), g("mad_gee_minus_mmm", "margin")),
             package = unlist(M[M$region == r, c("n_dp", "mad_gee", "mad_mmm", "diff", "margin")])) }))
same_classes <- setequal(M$region, H$region[grepl("^map\\|", H$region) & H$q == "mad_gee_minus_mmm"]) &&
  all(vapply(M$region, function(r) M$class[M$region == r] == H$class[H$region == r & H$q == "mad_gee_minus_mmm"], logical(1)))
cmp <- rbind(cmp, cmp2, cmp3)
cmp$equal <- eq(cmp$file, cmp$package)
other <- c(design_points = st$design_points == unique(T4$n_dp), all_strategies = setequal(both, T4$strategy),
           chosen = st$chosen == T4$strategy[T4$chosen], ar1_fg_points = same_points, map_classes = same_classes,
           md5_inputs = all(md5_now == MD5))
txt <- c(list.files("R", full.names = TRUE), list.files("man", full.names = TRUE), list.files("vignettes", "Rmd$", full.names = TRUE),
         "README.md", "DESCRIPTION", "inst/CITATION")
stale <- do.call(rbind, lapply(txt, function(f) { l <- readLines(f, warn = FALSE)
  k <- grep("1,896|1896|manuscript v4|slightly closer to 5%|recommended strategy, at every number of clusters and occasions,$|occasions, is generalized estimating|^ *is generalized estimating equations with", l)
  if (length(k)) data.frame(file = f, line = k, text = substr(l[k], 1, 120)) else NULL }))
write.csv(cmp, file.path(OUT, "c5_constants.csv"), row.names = FALSE)
sink(file.path(OUT, "c5_summary.txt"), split = TRUE)
cat("md5 of the input files:\n"); print(data.frame(file = names(MD5), plan = MD5, now = unname(md5_now)), row.names = FALSE)
cat("values compared:", nrow(cmp), "; equal:", sum(cmp$equal), "; per source:\n"); print(table(cmp$what, cmp$equal))
cat("strategies:", length(both), "of", nrow(T4), "; chosen in the file:", T4$strategy[T4$chosen], "; in the package:", st$chosen, "\n")
print(other); if (!all(cmp$equal)) print(cmp[!cmp$equal, ])
cat("lines with the old recommendation:", if (is.null(stale)) 0L else nrow(stale), "\n"); if (!is.null(stale)) print(stale)
cat("version in DESCRIPTION:", read.dcf("DESCRIPTION")[1, "Version"], "; NEWS opens with 0.3.1:", grepl("0.3.1", readLines("NEWS.md")[1]), "\n")
ok <- all(cmp$equal) && all(other) && is.null(stale) && read.dcf("DESCRIPTION")[1, "Version"] == "0.3.1" && grepl("0.3.1", readLines("NEWS.md")[1])
cat("GATE check 9:", if (ok) "PASS" else "FAIL", "\n"); sink()
