#' Tuning constants for the fits
#'
#' Most users never need this function: the defaults are those of the
#' simulation study behind the package.  It is passed to [fit_fewclust()]
#' as `control`.
#'
#' @param q Number of Gauss-Hermite quadrature points of the marginalized
#'   multilevel model.
#' @param check Whether the marginalized model applies its convergence
#'   check.  When `time` is given the model lets the latent standard
#'   deviation follow a linear trend over the occasions; a fit whose
#'   optimiser did not converge (`nlm()` code other than 1 or 2), or whose
#'   slope of the latent standard deviation is below `check_threshold`, is
#'   replaced by the fit with a constant latent standard deviation, but only
#'   if `nlm()` reports convergence of that fit (code 1 or 2) and its
#'   estimates and standard errors are finite.  Otherwise the marginalized
#'   model gives no result: [fit_fewclust()] falls back to GEE with
#'   `method = "auto"` and stops with an error with `method = "mmm"`.  The
#'   fit at the boundary is never reported.  When `time` is `NULL` the
#'   check likewise rejects a constant-SD fit that does not converge or has
#'   non-finite estimates or standard errors.  The model is recommended only
#'   with the check; switching it off is recorded in the notes of the fit.
#' @param check_coherence Whether the check also replaces a trend fit whose
#'   log-likelihood is below that of the fit with a constant latent standard
#'   deviation.  The constant fit is a sub-model of the trend fit, so a
#'   maximum cannot lie below it: the optimiser stopped at a false optimum.
#'   This condition is an addition of the package to the procedure of the
#'   study; `check_coherence = FALSE` gives that procedure exactly.
#' @param check_threshold Slope of the latent standard deviation (per
#'   occasion) below which the trend fit is considered at the boundary.
#' @param fg_bound Bound \eqn{b} of the Fay-Graubard correction, strictly
#'   inside (0, 1); 0.75 is the value recommended by Fay and Graubard (2001).
#'
#' @return A list of class `"fewclust_control"`.
#' @examples
#' fewclust_control()
#'
#' # the convergence check exactly as in the simulation study
#' fit_fewclust(fc_trial, y ~ time * arm,
#'   id = id, time = time,
#'   control = fewclust_control(check_coherence = FALSE)
#' )
#' @export
fewclust_control <- function(q = 20L, check = TRUE, check_coherence = TRUE,
                             check_threshold = 1e-3, fg_bound = 0.75) {
  if (!fc_is_whole(q) || q < 2) {
    fc_abort("{.arg q} must be a single whole number of at least 2.")
  }
  if (!fc_is_flag(check)) {
    fc_abort("{.arg check} must be TRUE or FALSE.")
  }
  if (!fc_is_flag(check_coherence)) {
    fc_abort("{.arg check_coherence} must be TRUE or FALSE.")
  }
  if (!fc_is_number(check_threshold) || check_threshold <= 0) {
    fc_abort("{.arg check_threshold} must be a single positive number.")
  }
  if (!fc_is_number(fg_bound) || fg_bound <= 0 || fg_bound >= 1) {
    fc_abort("{.arg fg_bound} must be a single number strictly inside (0, 1).")
  }
  structure(
    list(
      q = as.integer(q), check = check, check_coherence = check_coherence,
      check_threshold = check_threshold, fg_bound = fg_bound
    ),
    class = "fewclust_control"
  )
}

#' @export
print.fewclust_control <- function(x, ...) {
  cat("<fewclust_control>\n")
  cat(sprintf("  quadrature points : %d\n", x$q))
  cat(sprintf(
    "  convergence check : %s (threshold %g%s)\n",
    if (x$check) "on" else "off", x$check_threshold,
    if (x$check_coherence) ", with the coherence condition" else ""
  ))
  cat(sprintf("  Fay-Graubard bound: %g\n", x$fg_bound))
  invisible(x)
}
