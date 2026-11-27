# The engine follows the matrix notation of its references (X, Z, W, V,
# ...), so object names are not snake_case in this file.
# nolint start: object_name_linter.

# Pre-fit separation detection on the fixed-effects design.
#
# Two passes, both run BEFORE any clustered model is fitted:
#
#  1. per-column check (human-readable evidence): for a single column (with
#     intercept), detects strict range separation (complete) and boundary-tie
#     overlap (quasi-complete, reported at column level for discrete columns
#     with at most `max_tab_levels` distinct values so the message names the
#     boundary value).
#     A pure cell that is NOT extreme (e.g. all events at the top dose with
#     mixed outcomes below) does not make the MLE infinite and is not flagged.
#
#  2. exact linear program (Konis 2007): checks whether there exists a
#     non-zero separating vector beta in [-1, 1]^p with
#     (2 * y_i - 1) * (x_i' beta) >= 0 for all i and strictly positive for
#     at least one i.  Columns are scaled by max(abs(col)) so the box
#     constraint is invariant to covariate units (with or without an
#     intercept).  When separation is found, a second LP maximises the
#     minimum margin to distinguish complete from quasi-complete separation.

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

  lp_res <- fc_separation_lp(X, y)
  if (lp_res$separated) {
    quasi <- TRUE
    if (lp_res$complete) complete <- TRUE
    if (!length(evidence)) {
      evidence <- c(evidence, lp_res$evidence)
    }
  }

  status <- if (complete) "complete" else if (quasi) "quasi" else "none"
  list(status = status, evidence = evidence)
}

fc_separation_lp <- function(X, y) {
  none <- list(separated = FALSE, complete = FALSE, evidence = character(0))
  n <- nrow(X)
  p <- ncol(X)
  if (n == 0L || p == 0L) {
    return(none)
  }
  # Scale each non-zero column by max(abs) so [-1, 1]^p is unit-invariant
  col_max <- apply(abs(X), 2L, max)
  col_max[col_max == 0] <- 1
  Xs <- X / rep(col_max, each = n)
  Xbar <- Xs * (2 * y - 1)

  # Shift beta = z - 1 with 0 <= z <= 2 so lpSolve uses non-negative vars
  obj <- colSums(Xbar)
  rhs_sep <- rowSums(Xbar)
  mat <- rbind(Xbar, diag(p))
  dir <- c(rep(">=", n), rep("<=", p))
  rhs <- c(rhs_sep, rep(2, p))

  sol <- tryCatch(
    lpSolve::lp("max", obj, mat, dir, rhs),
    error = function(e) NULL
  )
  if (is.null(sol) || sol$status != 0L) {
    return(none)
  }
  beta <- sol$solution - 1
  margin <- drop(Xbar %*% beta)
  if (min(margin) < -1e-7 || max(margin) <= 1e-6 || sum(margin) <= 1e-6) {
    return(none)
  }

  is_complete <- all(margin > 1e-6)
  if (!is_complete) {
    # Maximize minimum margin t in [0, 1] subject to Xbar (z - 1) - t >= 0
    mat_c <- rbind(
      cbind(Xbar, rep(-1, n)),
      diag(p + 1L)
    )
    dir_c <- c(rep(">=", n), rep("<=", p + 1L))
    rhs_c <- c(rhs_sep, rep(2, p), 1)
    obj_c <- c(numeric(p), 1)
    sol_c <- tryCatch(
      lpSolve::lp("max", obj_c, mat_c, dir_c, rhs_c),
      error = function(e) NULL
    )
    if (!is.null(sol_c) && sol_c$status == 0L && sol_c$objval > 1e-6) {
      beta_c <- sol_c$solution[seq_len(p)] - 1
      if (min(drop(Xbar %*% beta_c)) > 1e-6) {
        is_complete <- TRUE
      }
    }
  }

  n_strict <- sum(margin > 1e-6)
  ev <- if (is_complete) {
    paste0(
      "a linear combination of the covariates separates the outcome classes ",
      "completely"
    )
  } else {
    sprintf(
      paste0(
        "a linear combination of the covariates separates %d of %d ",
        "observations (quasi-complete)"
      ),
      n_strict, n
    )
  }
  list(separated = TRUE, complete = is_complete, evidence = ev)
}
# nolint end
