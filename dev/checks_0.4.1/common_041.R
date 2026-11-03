# common_041.R -- helpers shared by the checks of fewclust 0.4.1 (gates A and B): the changed set predicted from a
# record of 0.4.0, and the rule of the study (D-90) applied to the fits recorded by the package.
cols_id <- c("term", "estimate", "std.error", "statistic", "df", "p.value")
# 0.4.0 result that D-90 changes: the marginalized model was reported with the refit not converged, or with the
# boundary fit kept because the refit was not finite
predicted_change <- function(x) {
  if (isTRUE(x$refused) || !identical(x$method, "mmm") || is.null(x$convergence)) return(FALSE)
  ck <- x$convergence
  (identical(ck$path, "constant fit used") && !ck$constant$code %in% c(1L, 2L)) ||
    identical(ck$path, "constant fit failed, trend fit kept")
}
finite_rec <- function(r) !is.null(r) && all(is.finite(r$estimate)) && all(is.finite(r$std.error))
# path that D-90 gives on the fits recorded in a check record (NA when the check was not requested)
rule_path <- function(ck) {
  if (is.null(ck) || !isTRUE(ck$requested)) return(NA_character_)
  tr <- ck$trend; co <- ck$constant
  below <- finite_rec(co) && tr$logLik < co$logLik - 1e-4
  bnd <- !tr$code %in% c(1L, 2L) || tr$slope_min < ck$threshold || (isTRUE(ck$coherence) && below)
  if (!bnd) return("trend fit kept")
  if (finite_rec(co) && co$code %in% c(1L, 2L)) "constant fit used" else "no result"
}
same_table <- function(a, b) identical(unname(as.list(a[, cols_id])), unname(as.list(b[, cols_id])))
same_fits <- function(a, b) identical(a$trend, b$trend) && identical(a$constant, b$constant)
