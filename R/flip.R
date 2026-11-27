#' Sign-flip score test for one coefficient
#'
#' Tests that one coefficient of the mean model is zero with the sign-flip
#' score test of Hemerik, Goeman and Finos (2020), treating the groups as
#' independent units.  It gives a p-value only: no estimate and no
#' interval.  It is useful with 2 or 3 occasions and few groups, where the
#' models of [fit_fewclust()] often give no usable estimate and have little
#' power; there the test keeps its level, but it has little power too.
#'
#' The model without the tested term is fitted as a logistic regression
#' with independent observations ([stats::glm.fit()], tolerance 1e-12, at
#' most 100 iterations).  The effective score contribution of each group to
#' the tested term is its score contribution minus the projection on the
#' other coefficients.  The p-value compares the absolute sum of the
#' contributions with its distribution under changes of sign: all
#' \eqn{2^N} sign assignments when there are at most `enum_max` groups,
#' otherwise `reps` random assignments plus the identity.  With other
#' coefficients estimated the test is asymptotically exact and robust to
#' differences in variance between groups.
#'
#' The test is not computable, and returns `NA` with a warning that gives
#' the reason, when the outcome is constant; when the model without the
#' tested term fails, does not converge or diverges (a standardised
#' coefficient of 1,000 or more in absolute value, as with separation);
#' when its information matrix is singular; or when the score contributions
#' are all zero or numerically zero relative to the scale of the tested
#' covariate (the tested term is not identified, for instance when one arm
#' has only zero outcomes).
#'
#' @param data A data frame (or tibble).
#' @param formula Two-sided formula for the mean model, e.g.
#'   `y ~ time * arm`.
#' @param id The column that identifies the groups: their contributions
#'   change sign together.  A bare column name or a string.
#' @param time Optional numeric column of the measurement occasions.  It
#'   only orders the rows within a group, which makes the result identical
#'   whatever the order of the rows of `data`.
#' @param term The coefficient to test: a column name of the model matrix
#'   (for example `"time:arm"`).  By default the last one, which for
#'   `y ~ time * arm` is the interaction.
#' @param reps Number of random sign assignments when there are more than
#'   `enum_max` groups.
#' @param enum_max Largest number of groups for which all \eqn{2^N} sign
#'   assignments are enumerated.
#' @param seed Optional seed.  When given, the random signs are drawn after
#'   `set.seed(seed)` and the random number stream of the session is left
#'   untouched; when `NULL`, the current stream is used.
#'
#' @return An object of class `"fewclust_flip"`, a list with `p.value`
#'   (`NA` when not computable), `statistic` (the sum of the effective
#'   score contributions), `term`, `n_groups`, `n_signs` (sign assignments
#'   in the reference distribution, the identity included), `enumerated`,
#'   `computable`, `reason` (`NA` when computable), `contributions` (one
#'   per group) and `call`.  [tidy()] turns it into a one-row tibble.
#'
#' @references
#' Hemerik J, Goeman JJ, Finos L (2020). Robust testing in generalized
#' linear models by sign flipping score contributions. *Journal of the Royal
#' Statistical Society: Series B* 82(3), 841-864.
#' \doi{10.1111/rssb.12369}
#'
#' @examples
#' flip_test(fc_trial, y ~ time * arm, id = id, time = time, seed = 1)
#'
#' # with 3 occasions and 12 groups
#' d <- simulate_clustered_binary(n_groups = 12, n_occasions = 3, seed = 4)
#' d |>
#'   flip_test(y ~ time * arm, id = id, time = time, seed = 1) |>
#'   tidy()
#' @seealso [fit_fewclust()]
#' @export
flip_test <- function(data, formula, id, time = NULL, term = NULL,
                      reps = 1999, enum_max = 12, seed = NULL) {
  call <- rlang::current_env()
  cl <- match.call()
  if (!fc_is_whole(reps) || reps < 1) {
    fc_abort("{.arg reps} must be a single whole number of at least 1.")
  }
  if (!fc_is_whole(enum_max) || enum_max < 1 || enum_max > 20) {
    fc_abort("{.arg enum_max} must be a single whole number between 1 and 20.")
  }
  if (!is.null(seed) && !(is.numeric(seed) && length(seed) == 1L)) {
    fc_abort("{.arg seed} must be NULL or a single number.")
  }
  cols <- fc_columns(data, formula, rlang::enquo(id), rlang::enquo(time),
    id_missing = missing(id), call = call
  )
  if (length(fc_findbars(formula))) {
    fc_abort(c(
      "{.arg formula} must not contain random-effects terms.",
      "i" = paste0(
        "The sign-flip test treats the groups as independent units; ",
        "give them with {.arg id}."
      )
    ))
  }
  # a design that is not of full rank is not an error here: when the tested
  # column is zero or collinear with the others the test is not computable
  # and says why
  parsed <- fc_parse(formula, data, cols$id, cols$time,
    full_rank = FALSE, call = call
  )
  x_mat <- parsed$X
  if (ncol(x_mat) < 2L) {
    fc_abort("There is no coefficient to test besides the intercept.")
  }
  if (is.null(term)) term <- colnames(x_mat)[ncol(x_mat)]
  valid <- is.character(term) && length(term) == 1L
  if (!valid || !term %in% colnames(x_mat)) {
    fc_abort("{.arg term} must be one of: {.val {colnames(x_mat)}}.")
  }
  ord <- order(parsed$id)
  if (!is.null(parsed$time)) {
    ord <- order(parsed$id, parsed$data[[parsed$time]])
  }
  group <- as.integer(parsed$id[ord])
  y <- as.numeric(parsed$y[ord])
  k <- match(term, colnames(x_mat))
  x <- unname(x_mat[ord, k])
  x0 <- unname(x_mat[ord, -k, drop = FALSE])
  offset <- parsed$offset
  if (!is.null(offset)) offset <- offset[ord]
  n_groups <- nlevels(parsed$id)

  if (!is.null(seed)) {
    old <- globalenv()$.Random.seed
    on.exit(fc_restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  enumerated <- n_groups <= enum_max
  # the sign matrix is drawn first: the random stream does not depend on
  # whether the test turns out computable
  signs <- if (enumerated) {
    fc_all_signs(n_groups)
  } else {
    matrix(sample(c(-1, 1), reps * n_groups, replace = TRUE), reps, n_groups)
  }
  # a constant outcome gives no usable result for any method
  res <- if (all(y == y[1L])) {
    list(
      p = NA_real_,
      why = paste0(
        "constant outcome (all ", y[1L], "): no usable ",
        "result, as for every strategy"
      ),
      stat = NA_real_, u = NULL
    )
  } else {
    fc_flip_core(group, y, x, x0, offset, signs, enumerated)
  }
  out <- structure(
    list(
      p.value = res$p, statistic = res$stat, term = term,
      n_groups = n_groups,
      n_signs = if (enumerated) nrow(signs) else nrow(signs) + 1L,
      enumerated = enumerated, computable = is.na(res$why),
      reason = res$why, contributions = res$u, call = cl
    ),
    class = "fewclust_flip"
  )
  if (!out$computable) {
    fc_warn(c(
      "The sign-flip test is not computable: {res$why}.",
      "i" = "The p-value is {.val {NA}}."
    ))
  }
  out
}

# Put back the random number stream of the session.
fc_restore_seed <- function(old) {
  if (is.null(old)) {
    rm(".Random.seed", envir = globalenv())
  } else {
    env <- globalenv()
    env[[".Random.seed"]] <- old
  }
}

# All 2^n sign vectors, first row all +1 (the identity).
fc_all_signs <- function(n) {
  m <- as.matrix(expand.grid(rep(list(c(1, -1)), n)))
  dimnames(m) <- NULL
  m
}

# The test at b0 = 0 with the non-computability rules of the study
# (morris-v4: R/s7d2/s7d2_ci.R::flip_p_b with diverge = TRUE).
fc_flip_core <- function(id, y, x, x0, offset, signs, enumerated) {
  out <- function(p = NA_real_, why = NA_character_, stat = NA_real_,
                  u = NULL) {
    list(p = p, why = why, stat = stat, u = u)
  }
  # Normalise nuisance columns by RMS so the reduced fit, its divergence
  # guard and its information matrix are invariant to covariate units.
  s0 <- apply(x0, 2L, function(v) sqrt(mean(v^2)))
  s0[!is.finite(s0) | s0 == 0] <- 1
  x0s <- x0 / rep(s0, each = nrow(x0))

  ctrl <- list(epsilon = 1e-12, maxit = 100, trace = FALSE)
  f <- tryCatch(
    suppressWarnings(
      glm.fit(x0s, y, family = binomial(), offset = offset, control = ctrl)
    ),
    error = function(e) NULL
  )
  ok <- !is.null(f) && isTRUE(f$converged) &&
    all(is.finite(f$coefficients)) && all(is.finite(f$fitted.values))
  if (!ok) {
    return(out(why = "reduced fit failed"))
  }
  if (max(abs(f$coefficients)) >= 1e3) {
    return(out(why = "reduced fit diverged"))
  }
  mu <- f$fitted.values
  e <- y - mu
  w <- mu * (1 - mu)
  info_00 <- crossprod(x0s, x0s * w)
  info_x0 <- colSums(x0s * (w * x))
  proj <- tryCatch(solve(info_00, info_x0), error = function(e) NULL)
  if (is.null(proj)) {
    return(out(why = "information singular"))
  }
  u_raw <- rowsum(x * e, id, reorder = TRUE)[, 1]
  u <- as.numeric(u_raw - rowsum(x0s * e, id, reorder = TRUE) %*% proj)
  if (!(sum(u^2) > 0)) {
    return(out(why = "zero score variance", stat = sum(u), u = u))
  }
  sx <- sqrt(mean(x^2))
  if (!is.finite(sx) || sx == 0) sx <- 1
  if (max(abs(u)) < 1e-6 * sx) {
    return(out(why = "score numerically zero", stat = sum(u), u = u))
  }
  t_obs <- abs(sum(u))
  t_flip <- abs(as.numeric(signs %*% u))
  p <- if (enumerated) {
    mean(t_flip >= t_obs * (1 - 1e-10))
  } else {
    (1 + sum(t_flip >= t_obs * (1 - 1e-10))) / (nrow(signs) + 1)
  }
  out(p = p, stat = sum(u), u = u)
}

#' @export
print.fewclust_flip <- function(x, digits = 4, ...) {
  cat("Sign-flip score test (p-value only; no estimate, no interval)\n\n")
  cat(sprintf("  H0        : coefficient of %s = 0\n", x$term))
  reference <- if (x$enumerated) {
    sprintf("all %s sign assignments", format(x$n_signs, big.mark = ","))
  } else {
    sprintf(
      "%s random sign assignments plus the identity",
      format(x$n_signs - 1L, big.mark = ",")
    )
  }
  cat(sprintf(
    "  groups    : %d; reference distribution: %s\n", x$n_groups,
    reference
  ))
  if (x$computable) {
    cat(sprintf(
      "  statistic : %s (sum of the effective score contributions)\n",
      format(x$statistic, digits = digits)
    ))
    cat(sprintf("  p-value   : %s\n", format.pval(x$p.value, digits = digits)))
  } else {
    cat(sprintf("  not computable: %s\n", x$reason))
  }
  invisible(x)
}

#' @rdname tidy.fewclust
#' @export
tidy.fewclust_flip <- function(x, ...) {
  tibble::tibble(
    term = x$term, statistic = x$statistic, p.value = x$p.value,
    n_groups = x$n_groups, n_signs = x$n_signs,
    method = "sign-flip score test"
  )
}
