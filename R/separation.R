# The engine follows the matrix notation of its references (X, Z, W, V,
# ...), so object names are not snake_case in this file.
# nolint start: object_name_linter.

# Pre-fit separation detection on the fixed-effects design.
#
# Two passes, both run BEFORE any clustered model is fitted:
#
#  1. per-column exact check: for a single column (with intercept), the MLE
#     is infinite if and only if the column separates the two outcome
#     classes -- strictly (no overlap of ranges: complete separation) or
#     with overlap at exactly one boundary value (quasi-complete).  The
#     strict check is exact for any column; the boundary-tie check is only
#     meaningful for discrete covariates and is therefore restricted to
#     columns with at most `max_tab_levels` (12) distinct values.
#     A pure cell that is NOT extreme (e.g. all events at the top dose with
#     mixed outcomes below) does not make the MLE infinite and is correctly
#     not flagged.
#
#  2. divergence probe: an unpenalised glm.fit run stopped at 25 and again
#     at 100 iterations; a linear predictor that keeps growing between the
#     two stops, together with fitted values at the boundary and large
#     standardised coefficients, marks separation that pass 1 cannot see
#     (a covariate combination).
#
# Probe-based detection of multi-column separation is conservative by
# design; the exact linear-programming test of Konis (2007) is named in the
# refusal message as the reference method.

detect_separation <- function(X, y, max_tab_levels = 12L,
                              has_intercept = "(Intercept)" %in%
                                colnames(X)) {
  if (length(unique(y)) < 2L) {
    return(list(status = "none", evidence = character(0)))
  }
  evidence <- character(0)
  complete <- FALSE
  quasi <- FALSE

  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    n_distinct <- length(unique(xj))
    if (n_distinct < 2L) next
    x1 <- xj[y == 1L]
    x0 <- xj[y == 0L]
    if (has_intercept) {
      # with an intercept, any separating threshold works
      sep_complete <- min(x1) > max(x0) || min(x0) > max(x1)
      sep_quasi <- min(x1) >= max(x0) || min(x0) >= max(x1)
      cut <- if (min(x1) >= max(x0)) max(x0) else max(x1)
    } else {
      # without an intercept the threshold is pinned at zero: disjoint
      # strictly-positive ranges are NOT separation (the MLE is finite)
      sep_complete <- (min(x1) > 0 && max(x0) < 0) ||
        (max(x1) < 0 && min(x0) > 0)
      sep_quasi <- (min(x1) >= 0 && max(x0) <= 0) ||
        (max(x1) <= 0 && min(x0) >= 0)
      cut <- 0
    }
    if (sep_complete) {
      complete <- TRUE
      evidence <- c(evidence, sprintf(paste0(
        "column '%s' separates the outcome classes completely (their value ",
        "ranges do not overlap)"
      ), colnames(X)[j]))
    } else if (n_distinct <= max_tab_levels && sep_quasi) {
      quasi <- TRUE
      evidence <- c(evidence, sprintf(paste0(
        "column '%s' separates the outcome classes except at the boundary ",
        "value %s (quasi-complete)"
      ), colnames(X)[j], format(cut)))
    }
  }

  probe <- fc_separation_probe(X, y, center = has_intercept)
  if (probe$separated) {
    quasi <- TRUE
    if (probe$all_boundary) complete <- TRUE
    evidence <- c(evidence, probe$evidence)
  }

  status <- if (complete) "complete" else if (quasi) "quasi" else "none"
  list(status = status, evidence = evidence)
}

fc_separation_probe <- function(X, y, center = TRUE) {
  none <- list(
    separated = FALSE, all_boundary = FALSE,
    evidence = character(0)
  )
  # standardise non-constant columns so the coefficient-size criterion has
  # a scale (X is full rank: fc_parse rejects rank-deficient designs);
  # centering is only model-preserving when an intercept is present
  Xs <- X
  varying <- apply(Xs, 2L, function(v) length(unique(v)) > 1L)
  if (any(varying)) {
    Xs[, varying] <- scale(Xs[, varying, drop = FALSE],
      center = center, scale = TRUE
    )
  }
  f1 <- tryCatch(
    suppressWarnings(glm.fit(Xs, y,
      family = binomial(),
      control = list(maxit = 25L)
    )),
    error = function(e) NULL
  )
  if (is.null(f1)) {
    return(none)
  }
  f2 <- tryCatch(
    suppressWarnings(glm.fit(Xs, y,
      family = binomial(),
      start = f1$coefficients,
      control = list(maxit = 100L)
    )),
    error = function(e) NULL
  )
  if (is.null(f2)) {
    return(none)
  }

  grew <- max(abs(f2$linear.predictors)) - max(abs(f1$linear.predictors)) > 2
  at_boundary <- f2$fitted.values < 1e-7 | f2$fitted.values > 1 - 1e-7
  big <- max(abs(f2$coefficients[varying]), 0) > 10
  if (any(at_boundary) && grew && big) {
    list(
      separated = TRUE,
      all_boundary = all(at_boundary),
      evidence = sprintf(
        paste0(
          "diverging fit: max |linear predictor| grew from %.1f to %.1f ",
          "between 25 and 100 iterations, with %d fitted value(s) at the ",
          "0/1 boundary"
        ),
        max(abs(f1$linear.predictors)), max(abs(f2$linear.predictors)),
        sum(at_boundary)
      )
    )
  } else {
    none
  }
}
# nolint end
