# c2_study_compare.R -- fewclust 0.4.1, gates A, B, C on the study sample and predictions 1-3
# (dev/aggiornamento_0.4.1.md). Reads the two files of c2_study_run.R (0.4.0 and 0.4.1) and the sample.
# Usage: Rscript --vanilla c2_study_compare.R <st_040.rds> <st_041.rds> <sample.csv> <md5_sources_before.csv> <out_dir>
args <- commandArgs(TRUE); OUT <- args[5]
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common_041.R"))
a <- readRDS(args[1]); b <- readRDS(args[2]); S <- fread(args[3]); M5 <- fread(args[4])
stopifnot(a$version == "0.4.0", b$version == "0.4.1", length(a$res) == nrow(S), length(b$res) == nrow(S))
eq <- function(x, y, tol) isTRUE(all.equal(x, y, tolerance = tol))
term_of <- function(x, what) if (!identical(x$status, "fitted")) NA_real_ else x$table[[what]][match("time:trt", x$table$term)]
ck_path <- function(x) if (is.null(x$convergence)) NA_character_ else x$convergence$path
rows <- vector("list", nrow(S))
for (k in seq_len(nrow(S))) {
  s <- S[k]; x <- a$res[[k]]; y <- b$res[[k]]
  if (!is.null(x$error) || !is.null(y$error) || !isTRUE(y$generated)) {
    rows[[k]] <- data.table(k = k, check_error = paste(x$error, y$error, if (!isTRUE(y$generated)) "not generated"))
    next
  }
  fx <- x$fits; fy <- y$fits
  # gate A: identity outside the predicted set, expected change inside, for mmm, auto and the engine
  as40 <- function(f) list(refused = f$status != "fitted", method = f$method, convergence = f$convergence)
  pred <- c(mmm = predicted_change(as40(fx$mmm)), auto = predicted_change(as40(fx$auto)),
            engine = identical(fx$engine$status, "fitted") && predicted_change(list(refused = FALSE, method = "mmm",
                                                                   convergence = fx$engine$convergence)))
  same_api <- function(p, q) identical(p$status, q$status) && (p$status != "fitted" ||
    (identical(p$method, q$method) && same_table(p$table, q$table)))
  same <- c(mmm = same_api(fx$mmm, fy$mmm), auto = same_api(fx$auto, fy$auto), gee = same_api(fx$gee, fy$gee),
            engine = identical(fx$engine$status, fy$engine$status) && identical(fx$engine$est, fy$engine$est) &&
              identical(fx$engine$se, fy$engine$se))
  changed_ok <- c(
    mmm = !pred[["mmm"]] || (fy$mmm$status == "fewclust_error" && identical(ck_path(fy$mmm), "no result") &&
                               grepl("did not converge \\(nlm code|no finite estimate", gsub("\\s+", " ", fy$mmm$message))),
    # auto: GEE identical to method = "gee"; where GEE fails too, the error of "neither model", naming the reason
    auto = !pred[["auto"]] || (same_api(fy$auto, fy$gee) && identical(ck_path(fy$auto), "no result") &&
      (if (fy$auto$status == "fitted") fy$auto$method == "gee" else
        grepl("did not converge \\(nlm code|no finite estimate", gsub("\\s+", " ", fy$auto$message)))),
    engine = !pred[["engine"]] || (fy$engine$status == "failure" && identical(ck_path(fy$engine), "no result")))
  gateA <- all(same[c("mmm", "auto", "engine")] | pred) && same[["gee"]] && all(changed_ok)
  rec_same <- if (!is.null(fx$engine$convergence) && !is.null(fy$engine$convergence))
    same_fits(fx$engine$convergence, fy$engine$convergence) else NA
  # gate B: the path of 0.4.1 is the rule of D-90 applied to the package's own recorded fits, and the outcome follows it
  rule_ok <- function(f, kind) {
    ck <- f$convergence; if (is.null(ck)) return(TRUE)
    rp <- rule_path(ck); if (is.na(rp)) return(TRUE)
    if (rp != ck$path) return(FALSE)
    if (rp != "no result") return(TRUE)
    switch(kind, mmm = f$status == "fewclust_error",
           auto = (f$status == "fitted" && f$method == "gee") || (f$status == "fewclust_error" && fy$gee$status != "fitted"),
           engine = f$status == "failure")
  }
  gateB <- rule_ok(fy$mmm, "mmm") && rule_ok(fy$auto, "auto") && rule_ok(fy$engine, "engine")
  ce <- fy$engine$convergence
  rows[[k]] <- data.table(k = k, check_error = NA_character_, data_identical = y$data_identical, n_obs = y$n_obs,
    sum_y = y$sum_y, gateA_ok = gateA, gateB_ok = gateB, recorded_fits_identical = rec_same,
    pred_change_mmm = pred[["mmm"]], pred_change_engine = pred[["engine"]],
    mmm_status = fy$mmm$status, mmm_path = ck_path(fy$mmm), mmm_message = if (is.null(fy$mmm$message)) NA_character_ else
      gsub("\\s+", " ", fy$mmm$message), mmm_est = term_of(fy$mmm, "estimate"), mmm_se = term_of(fy$mmm, "std.error"),
    auto_method = if (fy$auto$status == "fitted") fy$auto$method else NA_character_,
    engine_status = fy$engine$status, engine_path = if (is.null(ce)) NA_character_ else ce$path,
    engine_refit_status = if (is.null(ce)) NA_character_ else ce$refit_status,
    pkg_trend_code = if (is.null(ce$trend)) NA_integer_ else ce$trend$code,
    pkg_trend_slope = if (is.null(ce$trend)) NA_real_ else ce$trend$slope_min,
    pkg_const_code = if (is.null(ce$constant)) NA_integer_ else ce$constant$code,
    pkg_const_est = if (is.null(ce$constant)) NA_real_ else unname(ce$constant$estimate["time:trt"]),
    pkg_const_se = if (is.null(ce$constant)) NA_real_ else unname(ce$constant$std.error["time:trt"]),
    pkg_const_finite = finite_rec(ce$constant),
    engine_est = if (is.null(fy$engine$est)) NA_real_ else fy$engine$est,
    engine_se = if (is.null(fy$engine$se)) NA_real_ else fy$engine$se)
}
R <- cbind(S[, .(dgp, scen, N, ntime, label, sim_id, kind, why, A_code, A_slope_min, B_code, B_est, B_se, B_ok, used)],
           rbindlist(rows, fill = TRUE))
R[, api_no_result := mmm_status != "fitted"]
R[, no_result_by_check := mmm_status == "fewclust_error" & mmm_path %in% "no result"]
R[, design_refusal := mmm_status == "fewclust_error" & is.na(mmm_path)]
R[, engine_no_result := engine_status != "fitted"]
R[, agree_B := mmm_status == "fitted" & mmm_path %in% "constant fit used" &
      mapply(eq, mmm_est, B_est, 1e-6) & mapply(eq, mmm_se, B_se, 1e-5)]
R[, agree_B_engine := engine_status == "fitted" & engine_path %in% "constant fit used" &
      mapply(eq, engine_est, B_est, 1e-6) & mapply(eq, engine_se, B_se, 1e-5)]
fwrite(R, file.path(OUT, "c2_study_rows.csv"))
md5_now <- tools::md5sum(M5$file)
srcs_ok <- identical(unname(md5_now), M5$md5) && identical(a$md5_before, a$md5_after) &&
  identical(b$md5_before, b$md5_after) && identical(a$md5_after, b$md5_before)
gA <- R[, all(is.na(check_error)) && all(gateA_ok) && all(recorded_fits_identical, na.rm = TRUE)]
gB <- R[, all(gateB_ok)]
gC <- R[, all(is.na(check_error)) && all(data_identical %in% TRUE)] && srcs_ok
EXC <- R[(kind == "key" & !api_no_result) | (kind == "conv" & !agree_B)]
fwrite(EXC, file.path(OUT, "c2_exceptions.csv"))
sink(file.path(OUT, "c2_summary.txt"), split = TRUE)
options(width = 220)
cat("fewclust", a$version, "vs", b$version, "; datasets", nrow(R), "; seconds", a$seconds, "and", b$seconds, "\n")
cat("check errors:", R[!is.na(check_error), .N], "; data identical to the saved run (n_obs, sum_y):",
    R[data_identical %in% TRUE, .N], "of", nrow(R), "; study files unchanged (md5):", srcs_ok, "\n")
cat("gate A: identical to 0.4.0 outside the predicted set and changed as expected inside:", R[, sum(gateA_ok, na.rm = TRUE)],
    "of", nrow(R), "; predicted to change (mmm):", R[, sum(pred_change_mmm, na.rm = TRUE)], "; (engine):",
    R[, sum(pred_change_engine, na.rm = TRUE)], "; recorded fits identical:", R[, sum(recorded_fits_identical, na.rm = TRUE)],
    "of", R[, sum(!is.na(recorded_fits_identical))], "\n")
cat("gate B: rule of D-90 on the package's own fits:", R[, sum(gateB_ok, na.rm = TRUE)], "of", nrow(R), "\n")
cat("\nKeys (study: no result under D-90), method = 'mmm', check_coherence = FALSE:\n")
print(R[kind == "key", .(n = .N, no_result = sum(api_no_result), by_check = sum(no_result_by_check),
                         design_refusal = sum(design_refusal), result = sum(!api_no_result),
                         engine_no_result = sum(engine_no_result), auto_gee = sum(auto_method %in% "gee")),
        keyby = .(scen, why)])
cat("\nConverged refits (study used fit B with nlm code 1-2): agreement of time:trt with fit B (est 1e-6, SE 1e-5):\n")
print(R[kind == "conv", .(n = .N, agree = sum(agree_B), constant_used = sum(mmm_path %in% "constant fit used"),
                          no_result_by_check = sum(no_result_by_check), design_refusal = sum(design_refusal),
                          trend_kept = sum(mmm_path %in% "trend fit kept"), engine_agree = sum(agree_B_engine)),
        keyby = .(scen, dgp)])
cat("\nexceptions (written to c2_exceptions.csv):", nrow(EXC), "\n")
print(EXC[, .(kind, scen, dgp, N, ntime, label, sim_id, why, A_code, A_slope = signif(A_slope_min, 3), B_code,
              B_est = signif(B_est, 5), mmm_status, mmm_path, pkg_trend_code, pkg_trend_slope = signif(pkg_trend_slope, 3),
              pkg_const_code, pkg_const_finite, mmm_est = signif(mmm_est, 5), mmm_se = signif(mmm_se, 4),
              B_se = signif(B_se, 4))])
cat("\nGATE A (study sample):", if (gA) "PASS" else "FAIL", "\n")
cat("GATE B (study sample):", if (gB) "PASS" else "FAIL", "\n")
cat("GATE C:", if (gC) "PASS" else "FAIL", "\n")
sink()
