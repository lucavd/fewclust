# make_constants.R -- writes R/constants.R from the final tables of step S11 of the study (branch morris-v4 of
# lucavd/MMM, working copy; decisions D-89, separated datasets stay results, and D-90, the refit of the convergence
# check counts only if nlm converges). Nothing is typed by hand: every value is read from a table, and every value of
# R/constants.R carries, in a comment, its source file, the row it comes from and the id(s) of the registry of the
# study's numbers (tables/V5_S7f_numbers.csv) that cite the same cell, or "no registry id". Checked by
# dev/checks_0.4.2/c1_constants.R.
# Inputs (read only): ~/MMM/tables/V5_S7f_table4.csv (per-strategy table of the manuscript), V5_S7f_rule_dp.csv
#   (design points: grid, small region, clearly-above flags), V5_S7f_cp_S7d_dp_rates.csv (t reference),
#   V5_S7f_h2h_boot.csv (Figure 1a), V5_S7f_numbers.csv (registry).
# Usage (from the package root): Rscript --vanilla dev/make_constants.R
TB <- path.expand("~/MMM/tables")
files <- c(t4 = "V5_S7f_table4.csv", dp = "V5_S7f_rule_dp.csv", cp = "V5_S7f_cp_S7d_dp_rates.csv",
           hb = "V5_S7f_h2h_boot.csv", reg = "V5_S7f_numbers.csv")
T4 <- read.csv(file.path(TB, files[["t4"]]), stringsAsFactors = FALSE)
D <- read.csv(file.path(TB, files[["dp"]]), stringsAsFactors = FALSE)
K <- read.csv(file.path(TB, files[["cp"]]), stringsAsFactors = FALSE)
H <- read.csv(file.path(TB, files[["hb"]]), stringsAsFactors = FALSE)
REG <- read.csv(file.path(TB, files[["reg"]]), stringsAsFactors = FALSE, colClasses = "character")
md5 <- tools::md5sum(file.path(TB, files))
# the registry cites these files with their current md5 (when it cites them at all)
for (f in files[c("t4", "dp", "cp", "hb")]) {
  m <- unique(REG$md5[REG$file == file.path("tables", f)])
  stopifnot(length(m) <= 1L, !length(m) || m == unname(md5[file.path(TB, f)]))
}
# registry filter: a JSON object of string keys and values, e.g. {"strategy": "flip"}
parse_filter <- function(s) {
  kv <- regmatches(s, gregexpr('"[^"]*": *"[^"]*"', s))[[1]]
  stats::setNames(sub('^"[^"]*": *"([^"]*)"$', "\\1", kv), sub('^"([^"]*)".*$', "\\1", kv))
}
REG$keys <- lapply(REG$filter, parse_filter)
# ids of the registry citing the cell (file, keys, field); the registry value must equal the table value
reg_ids <- function(file, keys, field, value) {
  j <- which(REG$file == file.path("tables", file) & REG$field == field)
  j <- j[vapply(REG$keys[j], function(k) length(k) == length(keys) && all(names(keys) %in% names(k)) &&
                  all(k[names(keys)] == unlist(keys)), logical(1))]
  for (i in j) {
    rv <- REG$value[i]
    ok <- if (is.character(value)) identical(rv, value) else isTRUE(all.equal(as.numeric(rv), value, tolerance = 1e-12))
    if (!ok) stop("registry ", REG$id[i], " = ", rv, " but the table gives ", value)
  }
  if (length(j)) paste(REG$id[j], collapse = ", ") else "no registry id"
}
num <- function(x) if (is.character(x)) paste0("\"", x, "\"") else
  if (is.integer(x) || x == round(x)) paste0(format(as.integer(x), scientific = FALSE), "L") else sprintf("%.15g", x)
# one value per line, with its source and registry ids
col_code <- function(name, values, notes, indent = "    ") {
  body <- paste0(indent, "  ", vapply(values, num, ""), ifelse(seq_along(values) < length(values), ",", ""),
                 " # ", notes)
  c(paste0(indent, name, " = c("), body, paste0(indent, "),"))
}
df_code <- function(df, notes_of, rn = NULL, indent = "    ") {
  out <- character(0)
  for (n in names(df)) out <- c(out, col_code(n, df[[n]], notes_of(n), indent))
  if (!is.null(rn)) out <- c(out, paste0(indent, "row.names = c(", paste0("\"", rn, "\"", collapse = ", "), "),"))
  c(out, paste0(indent, "stringsAsFactors = FALSE"))
}
# 1. per-strategy table: every numeric column, row names = strategy
numcols <- names(T4)[vapply(T4, is.numeric, logical(1))]
t4 <- T4[, numcols]
t4_notes <- function(n) vapply(seq_len(nrow(T4)), function(i)
  sprintf("%s | %s | %s", T4$strategy[i], files[["t4"]], reg_ids(files[["t4"]], list(strategy = T4$strategy[i]), n, t4[[n]][i])), "")
# 2. design points where GEE AR(1) + Fay-Graubard + t is clearly above 7.5%
hi <- D[D$strategy == "gee_ar1_fg" & D$hi, c("dgp", "scen", "N", "ntime", "pu")]
hi <- merge(hi, K[K$strategy == "gee_ar1_fg" & K$ref == "t", c("dgp", "scen", "N", "ntime", "n", "nd", "kk", "rate")])
hi <- hi[order(hi$dgp), ]
stopifnot(nrow(hi) == T4$hi_all[T4$strategy == "gee_ar1_fg"])
hi_notes <- function(n) {
  src <- if (n %in% c("n", "nd", "kk", "rate")) files[["cp"]] else files[["dp"]]
  vapply(seq_len(nrow(hi)), function(i) {
    key <- if (src == files[["cp"]]) list(dgp = hi$dgp[i], scen = hi$scen[i], N = as.character(hi$N[i]),
      ntime = as.character(hi$ntime[i]), strategy = "gee_ar1_fg", ref = "t") else list(dgp = hi$dgp[i],
      scen = hi$scen[i], N = as.character(hi$N[i]), ntime = as.character(hi$ntime[i]), strategy = "gee_ar1_fg")
    sprintf("%s %s N=%d T=%d | %s | %s", hi$dgp[i], hi$scen[i], hi$N[i], hi$ntime[i], src, reg_ids(src, key, n, hi[[n]][i]))
  }, "")
}
# 3. Figure 1a map: mean absolute distance of the type I error from 5% (points), GEE AR(1) + FG minus MMM with the check
reg <- unique(H$region[grepl("^map\\|", H$region)])
cell <- function(r, q, col) { x <- H[H$region == r & H$q == q, col]; stopifnot(length(x) == 1L); x }
h2h <- do.call(rbind, lapply(reg, function(r) data.frame(region = r, n_dp = cell(r, "mad_gee", "n_dp"),
  mad_gee = cell(r, "mad_gee", "est"), mad_mmm = cell(r, "mad_mmm", "est"), diff = cell(r, "mad_gee_minus_mmm", "est"),
  margin = cell(r, "mad_gee_minus_mmm", "margin"), class = cell(r, "mad_gee_minus_mmm", "class"), stringsAsFactors = FALSE)))
h2h_src <- list(region = c("mad_gee", "region"), n_dp = c("mad_gee", "n_dp"), mad_gee = c("mad_gee", "est"),
                mad_mmm = c("mad_mmm", "est"), diff = c("mad_gee_minus_mmm", "est"), margin = c("mad_gee_minus_mmm", "margin"),
                class = c("mad_gee_minus_mmm", "class"))
h2h_notes <- function(n) vapply(seq_len(nrow(h2h)), function(i) {
  q <- h2h_src[[n]][1]; f <- h2h_src[[n]][2]
  sprintf("region %s, q %s, column %s | %s | %s", h2h$region[i], q, f, files[["hb"]],
          if (n == "region") "no registry id" else reg_ids(files[["hb"]], list(region = h2h$region[i], q = q), f, h2h[[n]][i]))
}, "")
# 4. scalars, all from the design points of tables/V5_S7f_rule_dp.csv
dp <- unique(D[, c("dgp", "scen", "N", "ntime", "small")])
stopifnot(!anyDuplicated(dp[, c("dgp", "scen", "N", "ntime")]))
Ns <- sort(unique(dp$N)); stopifnot(all(D$small == (D$ntime <= max(D$ntime[D$small]) & D$N <= max(D$N[D$small]))))
design_points <- unique(T4$n_dp); stopifnot(length(design_points) == 1L, design_points == nrow(dp))
stopifnot(sum(dp$small) == unique(T4$n_small))
chosen <- T4$strategy[T4$chosen]; stopifnot(length(chosen) == 1L)
dpn <- sprintf("design points of %s", files[["dp"]])
sc <- list(
  min_groups = list(Ns[1], paste("smallest simulated number of clusters: min(N),", dpn)),
  next_groups = list(Ns[2], paste("next simulated number of clusters: second value of N,", dpn)),
  max_groups = list(max(Ns), paste("largest: max(N),", dpn)),
  max_occasions = list(max(dp$ntime), paste("longest simulated series: max(ntime),", dpn)),
  short_series = list(max(dp$ntime[dp$small]), paste("occasions at or below which ... (max(ntime) where small),", dpn)),
  few_groups = list(max(dp$N[dp$small]), paste("... with this many clusters or fewer (max(N) where small): the small region of pu_small,", dpn)),
  design_points = list(design_points, sprintf("n_dp of %s (= rows of %s)", files[["t4"]], dpn)),
  mechanisms = list(length(unique(dp$dgp)), paste("number of data-generating mechanisms: distinct dgp,", dpn)),
  chosen = list(chosen, sprintf("strategy with chosen = TRUE in %s", files[["t4"]])))
sc_code <- vapply(names(sc), function(n) sprintf("  %s = %s, # %s; %s", n, num(sc[[n]][[1]]), sc[[n]][[2]],
  if (n == "chosen") reg_ids(files[["t4"]], list(strategy = chosen), "chosen", "TRUE") else
    if (n == "design_points") reg_ids(files[["t4"]], list(strategy = chosen), "n_dp", design_points) else "no registry id"), "")
out <- c(
"# Numbers of the simulation study behind the package (Vedovelli et al.,",
"# manuscript), from the final tables of step S11 of the study: decisions",
"# D-89 (separated datasets stay results) and D-90 (the refit of the",
"# convergence check of the marginalized model counts only if nlm",
"# converges), and D-82 (a dataset whose outcome is constant gives no usable",
"# result for any method).  The vignette \"What the simulation study found\"",
"# reads them from here.",
"#",
"# GENERATED by dev/make_constants.R from the tables of branch morris-v4 of",
"# lucavd/MMM (working copy of step S11; md5 below) and checked against them",
"# by dev/checks_0.4.2/c1_constants.R; do not edit by hand.  Next to every",
"# value: the row it comes from, the source file, and the id(s) of the",
"# registry of the study's numbers (tables/V5_S7f_numbers.csv) that cite the",
"# same cell, or \"no registry id\".",
sprintf("#   %s  %s", unname(md5), basename(names(md5))),
"#",
"# `table4`: the per-strategy table of the manuscript (tables/V5_S7f_table4.csv),",
"# one row per strategy, with a t reference on N - p degrees of freedom (the",
"# likelihood-ratio and sign-flip tests with their own references):",
"#   hi_all   design points clearly above a 7.5% type I error (of n_dp)",
"#   lo_all   design points clearly below 2.5%",
"#   lo_N10   of which with 10 clusters (of n_N10)",
"#   pu       mean share of datasets with a usable result (%)",
"#   pu_small the same with 2-3 occasions and at most 25 clusters (%)",
"# Row names follow the study: gee_<working correlation>_<variance>,",
"# mmm_cc (marginalized model with the convergence check), glmer (Wald test",
"# of the mixed model), glmer_lrt, flip (sign-flip score test), shrink.",
"# `ar1_fg_hi`: the design points where GEE AR(1) + Fay-Graubard + t is",
"# clearly above 7.5% (tables/V5_S7f_rule_dp.csv), with the usable datasets",
"# (n of nd), the rejections under H0 (kk) and the rate (%) of",
"# tables/V5_S7f_cp_S7d_dp_rates.csv; pu is the usable share (%).",
"# `h2h`: Figure 1a (tables/V5_S7f_h2h_boot.csv, regions map|<clusters>|",
"# <occasions>): mean absolute distance of the type I error from 5% (points)",
"# of GEE AR(1) + Fay-Graubard and of the marginalized model, their",
"# difference, the margin, and the class of the block bootstrap (GEE or MMM:",
"# that strategy is closer to 5% by more than the margin).",
".fc_study <- list(",
unname(sc_code),
"  table4 = data.frame(", df_code(t4, t4_notes, T4$strategy), "  ),",
"  ar1_fg_hi = data.frame(", df_code(hi, hi_notes), "  ),",
"  h2h = data.frame(", df_code(h2h, h2h_notes), "  )",
")")
writeLines(out, "R/constants.R")
cat("R/constants.R written:", nrow(t4), "strategies,", length(numcols), "columns;", nrow(hi), "points;", nrow(h2h),
    "map cells; chosen", chosen, "; registry ids attached:", sum(!grepl("no registry id$", out[grepl(" # ", out)])), "of",
    sum(grepl(" # ", out)), "values\n")
