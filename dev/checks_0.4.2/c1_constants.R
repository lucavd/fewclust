# c1_constants.R -- fewclust 0.4.2, gates 1 and 4 (dev/aggiornamento_0.4.2.md). Reads .fc_study from the installed
# library and from R/constants.R of the previous version, re-reads the source tables independently of
# dev/make_constants.R, and checks: every value = its table cell (integers and text identical, reals all.equal 1e-12);
# every registry id written next to a value cites a registry value equal to it; md5 of the sources = header = registry;
# same structure as the previous version. Writes the list of changed constants (old -> new, registry ids).
# Usage (from the package root): Rscript --vanilla dev/checks_0.4.2/c1_constants.R <lib> <old_constants.R> <out_dir>
args <- commandArgs(TRUE); LIB <- args[1]; OLD <- args[2]; OUT <- args[3]
.libPaths(c(LIB, "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
st <- get(".fc_study", envir = asNamespace(loadNamespace("fewclust", lib.loc = LIB)))
old <- local({ source(OLD, local = TRUE); .fc_study })
TB <- path.expand("~/MMM/tables")
T4 <- read.csv(file.path(TB, "V5_S7f_table4.csv"), stringsAsFactors = FALSE)
D <- read.csv(file.path(TB, "V5_S7f_rule_dp.csv"), stringsAsFactors = FALSE)
K <- read.csv(file.path(TB, "V5_S7f_cp_S7d_dp_rates.csv"), stringsAsFactors = FALSE)
H <- read.csv(file.path(TB, "V5_S7f_h2h_boot.csv"), stringsAsFactors = FALSE)
REG <- read.csv(file.path(TB, "V5_S7f_numbers.csv"), stringsAsFactors = FALSE, colClasses = "character")
same <- function(a, b) if (is.character(a) || is.logical(a)) identical(as.character(a), as.character(b)) else
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-12))
bad <- character(0); note <- function(x) bad <<- c(bad, x)
# 1. table4
for (s in T4$strategy) for (n in names(st$table4)) if (!same(st$table4[s, n], T4[T4$strategy == s, n])) note(paste("table4", s, n))
if (!identical(rownames(st$table4), T4$strategy)) note("table4 row names")
# 2. ar1_fg_hi, recomputed
dp_hi <- D[D$strategy == "gee_ar1_fg" & D$hi, ]
if (nrow(dp_hi) != nrow(st$ar1_fg_hi)) note("ar1_fg_hi rows")
for (i in seq_len(nrow(st$ar1_fg_hi))) {
  r <- st$ar1_fg_hi[i, ]
  d <- dp_hi[dp_hi$dgp == r$dgp & dp_hi$scen == r$scen & dp_hi$N == r$N & dp_hi$ntime == r$ntime, ]
  k <- K[K$strategy == "gee_ar1_fg" & K$ref == "t" & K$dgp == r$dgp & K$scen == r$scen & K$N == r$N & K$ntime == r$ntime, ]
  if (nrow(d) != 1L || nrow(k) != 1L || !same(r$pu, d$pu) || !all(mapply(same, r[c("n", "nd", "kk", "rate")], k[c("n", "nd", "kk", "rate")])))
    note(paste("ar1_fg_hi", i))
}
# 3. h2h
for (i in seq_len(nrow(st$h2h))) {
  r <- st$h2h[i, ]; g <- function(q, col) H[H$region == r$region & H$q == q, col]
  ok <- same(r$n_dp, g("mad_gee", "n_dp")) && same(r$mad_gee, g("mad_gee", "est")) && same(r$mad_mmm, g("mad_mmm", "est")) &&
    same(r$diff, g("mad_gee_minus_mmm", "est")) && same(r$margin, g("mad_gee_minus_mmm", "margin")) &&
    same(r$class, g("mad_gee_minus_mmm", "class"))
  if (!ok) note(paste("h2h", r$region))
}
if (!identical(st$h2h$region, unique(H$region[grepl("^map\\|", H$region)]))) note("h2h regions")
# 4. scalars
dp <- unique(D[, c("dgp", "scen", "N", "ntime", "small")]); Ns <- sort(unique(dp$N))
sc <- list(min_groups = Ns[1], next_groups = Ns[2], max_groups = max(Ns), max_occasions = max(dp$ntime),
           short_series = max(dp$ntime[dp$small]), few_groups = max(dp$N[dp$small]), design_points = nrow(dp),
           mechanisms = length(unique(dp$dgp)), chosen = T4$strategy[T4$chosen])
for (n in names(sc)) if (!identical(as.vector(st[[n]]), as.vector(if (is.character(sc[[n]])) sc[[n]] else as.integer(sc[[n]])))) note(paste("scalar", n))
# 5. registry ids written in R/constants.R
src <- readLines("R/constants.R")
ids <- unlist(lapply(regmatches(src, regexpr("f:[^ ,]+(, f:[^ ,]+)*$", src)), function(x) strsplit(x, ", ")[[1]]))
vals <- unlist(lapply(grep("f:[^ ,]+(, f:[^ ,]+)*$", src, value = TRUE), function(l) {
  v <- sub('^ *("[^"]*"|[^ ,]+),? # .*$', "\\1", sub("^ *[a-z_0-9]+ = ", "", l)); n <- length(strsplit(sub(".*\\| ", "", sub("; ", " | ", l)), ", ")[[1]])
  rep(v, n) }))
stopifnot(length(ids) == length(vals))
reg_ok <- mapply(function(id, v) {
  rv <- REG$value[REG$id == id]; if (length(rv) != 1L) return(FALSE)
  v <- sub("L$", "", gsub('"', "", v))
  if (v %in% c("TRUE", "FALSE") || is.na(suppressWarnings(as.numeric(v)))) identical(rv, v) else isTRUE(all.equal(as.numeric(rv), as.numeric(v), tolerance = 1e-12))
}, ids, vals)
if (!all(reg_ok)) note(paste("registry", ids[!reg_ok]))
# 6. md5: sources = header = registry
files <- c("V5_S7f_table4.csv", "V5_S7f_rule_dp.csv", "V5_S7f_cp_S7d_dp_rates.csv", "V5_S7f_h2h_boot.csv", "V5_S7f_numbers.csv")
m <- unname(tools::md5sum(file.path(TB, files)))
hdr <- vapply(files, function(f) sub("^#   ([0-9a-f]{32})  .*$", "\\1", grep(paste0("  ", f, "$"), src, value = TRUE)[1]), "")
regm <- vapply(files[1:4], function(f) { x <- unique(REG$md5[REG$file == file.path("tables", f)]); if (length(x)) x else NA_character_ }, "")
md5_ok <- identical(unname(hdr), m) && all(is.na(regm) | regm == m[1:4])
if (!md5_ok) note("md5")
# 7. structure as in the previous version
str_ok <- identical(names(st), names(old)) && identical(names(st$table4), names(old$table4)) &&
  identical(rownames(st$table4), rownames(old$table4)) && identical(names(st$ar1_fg_hi), names(old$ar1_fg_hi)) &&
  identical(names(st$h2h), names(old$h2h)) && identical(st$h2h$region, old$h2h$region)
if (!str_ok) note("structure")
# 8. changed constants, with the registry ids of the new file
idline <- function(element, column, row) {
  start <- grep(paste0("^  ", element, " = "), src)
  if (is.na(column)) return(sub("^.*; ", "", src[start]))
  b <- start + grep(paste0("^    ", column, " = c\\($"), src[start:length(src)])[1] - 1L
  sub("^.*\\| ", "", src[b + row])
}
rows <- list()
for (n in names(st)) {
  if (!is.data.frame(st[[n]])) {
    if (!identical(st[[n]], old[[n]])) rows[[length(rows) + 1L]] <- data.frame(element = n, row = NA, column = NA,
      old = format(old[[n]]), new = format(st[[n]]), registry = idline(n, NA, NA))
    next
  }
  a <- old[[n]]; b <- st[[n]]
  key <- if (n == "table4") rownames(b) else if (n == "h2h") b$region else paste(b$dgp, b$scen, b$N, b$ntime)
  keya <- if (n == "table4") rownames(a) else if (n == "h2h") a$region else paste(a$dgp, a$scen, a$N, a$ntime)
  for (i in seq_len(nrow(b))) for (cc in names(b)) {
    j <- match(key[i], keya); ov <- if (is.na(j)) NA else a[j, cc]
    if (is.na(j) || !identical(ov, b[i, cc])) rows[[length(rows) + 1L]] <- data.frame(element = n, row = key[i], column = cc,
      old = if (is.na(j)) "(absent)" else format(ov, digits = 15), new = format(b[i, cc], digits = 15), registry = idline(n, cc, i))
  }
  for (k in setdiff(keya, key)) rows[[length(rows) + 1L]] <- data.frame(element = n, row = k, column = "(row)", old = "present", new = "(absent)", registry = "")
}
DF <- do.call(rbind, rows)
if (is.null(DF)) DF <- data.frame(element = character(0), row = character(0), column = character(0), old = character(0), new = character(0), registry = character(0))
write.csv(DF, file.path(OUT, "constants_diff.csv"), row.names = FALSE)
sink(file.path(OUT, "c1_summary.txt"), split = TRUE)
cat("fewclust", as.character(packageVersion("fewclust", lib.loc = LIB)), "; values checked against the tables: table4",
    length(as.matrix(st$table4)), ", ar1_fg_hi", length(as.matrix(st$ar1_fg_hi)), ", h2h", length(as.matrix(st$h2h)),
    ", scalars", length(sc), "\n")
cat("registry ids written next to values:", length(ids), "; equal to the registry value:", sum(reg_ok), "\n")
cat("md5 sources = header = registry:", md5_ok, "; structure as before:", str_ok, "\n")
cat("constants changed:", nrow(DF), "(", paste(names(table(DF$element)), table(DF$element), collapse = "; "), ")\n")
if (length(bad)) { cat("discrepancies:\n"); print(bad) }
cat("GATE 1:", if (!length(bad)) "PASS" else "FAIL", "\n")
sink()
