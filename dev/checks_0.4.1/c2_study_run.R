# c2_study_run.R -- fewclust 0.4.1, check C: regenerate the sampled study datasets (out/c2_sample.csv) with the study's
# generators and seeds (R/s7d2/s7d2_common.R of branch morris-v4, read only, as dev/checks_0.3.1/c9_mmm_repro.R) and fit
# them with the fewclust of <lib>: method = "mmm" and "auto" with fewclust_control(check_coherence = FALSE) (the rule
# of the study), and method = "gee" (AR(1), Fay-Graubard). Mean model of the study: binY ~ time * trt (M5: + dose).
# Data identity against n_obs and sum_y of results_v4/S7d2/rows. At most 5 forked workers (6 R processes).
# Usage (from the package root): Rscript --vanilla dev/checks_0.4.1/c2_study_run.R <lib> <sample.csv> <out.rds> [n]
args <- commandArgs(TRUE); FC_LIB <- normalizePath(args[1]); SMP <- args[2]; OUTF <- args[3]
NMAX <- if (length(args) >= 4L) as.integer(args[4]) else Inf
PROJ <- "/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
.libPaths(c(FC_LIB, PROJ))
source(path.expand("~/MMM/R/s7d2/s7d2_common.R"))   # generators, seeds, tasks (read only)
# the study code may attach the fewclust of ~/R/s7d_lib: detach it and load the one of <lib>
# the study code also resets .libPaths() and defines globals (LIB among them): use FC_LIB and set the paths again
.libPaths(c(FC_LIB, PROJ))
if ("package:fewclust" %in% search()) detach("package:fewclust", unload = TRUE)
if ("fewclust" %in% loadedNamespaces()) unloadNamespace("fewclust")
invisible(loadNamespace("fewclust", lib.loc = FC_LIB))
stopifnot(dirname(getNamespaceInfo("fewclust", "path")) == FC_LIB)
VER <- as.character(getNamespaceVersion("fewclust"))
suppressMessages(library(data.table)); setDTthreads(1L)
S <- fread(SMP); if (is.finite(NMAX)) S <- head(S, NMAX)
TA <- s7d2_tasks()
S[, chunk := as.integer((sim_id - 1L) %/% 100L + 1L)]
S <- merge(S, TA[, .(dgp, scen, N, ntime, label, chunk, model, file)], by = c("dgp", "scen", "N", "ntime", "label", "chunk"),
           all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(S$file), all(file.exists(unique(S$file))))
md5_before <- tools::md5sum(unique(S$file))
fc <- asNamespace("fewclust")
paper <- fc$fewclust_control(check_coherence = FALSE)
rec_of <- function(r) if (is.null(r)) NULL else r[c("estimate", "std.error", "code", "slope_min", "logLik")]
one_fit <- function(d, f, method) {
  r <- tryCatch(suppressWarnings(suppressMessages(fc$fit_fewclust(d, f, id = "id", time = "time", method = method,
         control = paper))), fewclust_error = function(e) e, error = function(e) structure(list(other = conditionMessage(e)), class = "other"))
  if (inherits(r, "other")) return(list(status = paste("other error:", r$other)))
  if (inherits(r, "error")) return(list(status = "fewclust_error", message = conditionMessage(r), convergence = r$convergence))
  list(status = "fitted", method = r$method, fallback = r$fallback, table = as.data.frame(generics::tidy(r)),
       convergence = r$convergence, code = if (r$method == "mmm") r$engine$code else NA_integer_)
}
# the engine alone (fc_parse + fit_mmm, as tests/testthat/test-mmm-check.R): the convergence check also on datasets
# that fit_fewclust() refuses for the design (separation, outcome varying in one group), where the study fitted anyway
engine_fit <- function(d, f) {
  r <- tryCatch({
    parsed <- fc$fc_parse(f, d, quote(id), quote(time))
    fc$fit_mmm(parsed, paper, "t", 0.05, fc$fc_diagnose(parsed)$df_t)
  }, error = function(e) structure(list(reason = paste("error:", conditionMessage(e))), class = "other"))
  if (inherits(r, "other")) return(list(status = r$reason))
  if (inherits(r, "fc_failure")) return(list(status = "failure", reason = r$reason, convergence = r$check))
  tab <- r$coefficients; j <- match("time:trt", tab$term)
  list(status = "fitted", est = tab$estimate[j], se = tab$std.error[j], convergence = r$model$check, code = r$model$code)
}
one <- function(k) {
  t <- S[k]; i <- t$sim_id
  gen <- s7d2_gen(t)
  saved <- readRDS(t$file)[sim_id == i]
  deg <- design_degenerate(t$model, t$N, t$ntime, i)
  set_rep_seed(t$model, i)
  d <- tryCatch(as.data.frame(gen(t$N, t$ntime, i, t$label)), error = function(e) e)
  base <- list(k = k, dgp = t$dgp, scen = t$scen, N = t$N, ntime = t$ntime, label = t$label, sim_id = i, kind = t$kind,
               why = t$why, design_degenerate = deg)
  if (inherits(d, "error")) return(c(base, generated = FALSE, gen_error = conditionMessage(d)))
  d <- d[order(d$id, d$time), ]
  f <- if (t$scen == "M5") binY ~ time * trt + dose else binY ~ time * trt
  c(base, generated = TRUE, data_identical = nrow(saved) == 1L && nrow(d) == saved$n_obs && sum(d$binY) == saved$sum_y,
    n_obs = nrow(d), sum_y = sum(d$binY),
    fits = list(list(mmm = one_fit(d, f, "mmm"), auto = one_fit(d, f, "auto"), gee = one_fit(d, f, "gee"),
                     engine = engine_fit(d, f))))
}
one_safe <- function(k) tryCatch(one(k), error = function(e) list(k = k, error = conditionMessage(e)))
t0 <- proc.time()[[3]]
R <- parallel::mclapply(seq_len(nrow(S)), one_safe, mc.cores = 5L, mc.preschedule = FALSE)
secs <- round(proc.time()[[3]] - t0)
md5_after <- tools::md5sum(unique(S$file))
saveRDS(list(version = VER, lib = FC_LIB, sample = S, res = R, seconds = secs, md5_before = md5_before,
             md5_after = md5_after), OUTF)
cat("fewclust", VER, ": datasets", length(R), "; errors in the check:", sum(vapply(R, function(x) !is.null(x$error), logical(1))),
    "; seconds", secs, "\n")
