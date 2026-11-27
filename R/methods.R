#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment

fc_corstr_label <- function(corstr) {
  c(
    ar1 = "AR(1)", exchangeable = "exchangeable",
    independence = "independence"
  )[[corstr]]
}

fc_correction_label <- function(correction) {
  c(
    fg = "Fay-Graubard correction", md = "Mancl-DeRouen correction",
    kc = "Kauermann-Carroll correction", none = "no small-sample correction"
  )[[correction]]
}

fc_method_label <- function(x) {
  if (x$method == "mmm") {
    return("Marginalized multilevel model")
  }
  sprintf(
    "GEE, %s working correlation, %s", fc_corstr_label(x$corstr),
    fc_correction_label(x$correction)
  )
}

fc_check_label <- function(m) fc_check_text(m$check, m$code)

# Text of the record of the convergence check.  `code` is the nlm code of
# the reported fit (NA when the marginalized model gave no result).
fc_check_text <- function(ck, code = NA_integer_) {
  if (!ck$requested) {
    if (is.null(ck$trend)) {
      return("latent SD constant (no `time`)")
    }
    return("latent SD with time trend; convergence check switched off")
  }
  if (identical(ck$path, "no result")) {
    problem <- if (identical(ck$refit_status, "not converged")) {
      sprintf("did not converge (nlm code %d)", ck$constant$code)
    } else {
      "gave no finite estimate"
    }
    return(paste0(
      "no result: convergence check intervened (trend fit at the ",
      "boundary) and the refit without the trend ", problem
    ))
  }
  outcome <- if (!ck$boundary) {
    "passed"
  } else {
    "intervened (trend fit at the boundary)"
  }
  code_text <- ""
  if (!is.na(code) && !code %in% c(1L, 2L)) {
    code_text <- sprintf("; not converged (nlm code %d)", code)
  }
  sprintf(
    "latent SD %s; convergence check %s%s", ck$latent_sd, outcome,
    code_text
  )
}

fc_header <- function(x) {
  cat(fc_method_label(x), "\n", sep = "")
  if (!is.null(x$fallback)) {
    cat(
      "  (the marginalized multilevel model gave no usable estimate: ",
      x$fallback, ")\n",
      sep = ""
    )
  }
  cat(sprintf("Formula: %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf(
    "Data:    %d observations in %d groups (%s), %g occasions%s\n",
    x$n_obs, x$n_groups, x$id, x$n_occasions,
    if (is.null(x$time)) " (no `time`)" else sprintf(" (%s)", x$time)
  ))
  cat(sprintf(
    "Tests:   t reference with %d df (groups - coefficients)\n", x$df
  ))
}

fc_print_notes <- function(notes) {
  if (!length(notes)) {
    return(invisible())
  }
  cat("\nNotes:\n")
  for (n in notes) {
    cat(strwrap(n, indent = 0, exdent = 4, initial = "  * "), sep = "\n")
  }
}

#' @export
print.fewclust <- function(x, digits = max(3L, getOption("digits") - 3L),
                           ...) {
  fc_header(x)
  cat("\n")
  tab <- as.data.frame(tidy(x, conf.int = TRUE))
  rownames(tab) <- tab$term
  tab <- tab[setdiff(names(tab), c("term", "df"))]
  print(format(tab, digits = digits), quote = FALSE)
  fc_print_notes(x$notes)
  invisible(x)
}

#' Summarise a fewclust fit
#'
#' @param object,x A `"fewclust"` fit, or its summary.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#' @return `summary()` returns an object of class `"summary.fewclust"`,
#'   printed with the coefficient table, the confidence intervals, details
#'   of the model (convergence check of the marginalized model, working
#'   correlation of GEE) and the notes.
#' @examples
#' fit <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
#' summary(fit)
#' @export
summary.fewclust <- function(object, ...) {
  structure(list(fit = object, table = tidy(object, conf.int = TRUE)),
    class = "summary.fewclust"
  )
}

#' @rdname summary.fewclust
#' @export
print.summary.fewclust <- function(x,
                                   digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  fit <- x$fit
  tab <- as.data.frame(x$table)
  fc_header(fit)
  cat("\nCoefficients:\n")
  cm <- cbind(
    Estimate = tab$estimate, `Std. Error` = tab$std.error,
    `t value` = tab$statistic, `Pr(>|t|)` = tab$p.value
  )
  rownames(cm) <- tab$term
  stats::printCoefmat(cm,
    digits = digits, signif.stars = FALSE,
    na.print = "NA"
  )
  cat(sprintf("\n%s%% confidence intervals:\n", format(100 * fit$conf_level)))
  ci <- cbind(lower = tab$conf.low, upper = tab$conf.high)
  rownames(ci) <- tab$term
  print(signif(ci, digits))
  cat("\nModel:\n")
  if (fit$method == "mmm") {
    cat(sprintf("  %s\n", fc_check_label(fit$engine)))
    cat(sprintf(
      "  log-likelihood %.3f; nlm code %d after %d iterations\n",
      fit$engine$logLik, fit$engine$code, fit$engine$iterations
    ))
  } else {
    a <- fit$engine$geese$alpha
    alpha <- ""
    if (length(a)) {
      alpha <- paste0(", alpha = ", toString(format(a, digits = digits)))
    }
    cat(sprintf("  working correlation %s%s\n", fit$corstr, alpha))
    if (!is.null(fit$convergence)) {
      cat(sprintf(
        "  marginalized multilevel model: %s\n",
        fc_check_text(fit$convergence)
      ))
    }
  }
  cat(sprintf(
    "  outcome: event \"%s\", prevalence %.3f\n", fit$event,
    fit$design$prevalence
  ))
  fc_print_notes(fit$notes)
  invisible(x)
}

#' @export
print.fewclust_mmm <- function(x, digits = 4, ...) {
  cat("<fewclust_mmm> marginalized multilevel model\n")
  cat(sprintf(
    "  log-likelihood %.3f; nlm code %d after %d iterations\n",
    x$logLik, x$code, x$iterations
  ))
  cat("  ", fc_check_label(x), "\n\n", sep = "")
  print(round(cbind(estimate = x$theta, std.error = fc_se(x$vcov)), digits))
  invisible(x)
}

#' @export
coef.fewclust <- function(object, ...) object$coefficients

#' @export
vcov.fewclust <- function(object, ...) object$vcov

#' @export
nobs.fewclust <- function(object, ...) object$n_obs

#' @export
formula.fewclust <- function(x, ...) x$formula

# Wald intervals on the t reference; NA for coefficients without inference.
fc_ci <- function(x, level) {
  tab <- x$table
  crit <- stats::qt(1 - (1 - level) / 2, df = x$df)
  ok <- !is.na(tab$statistic)
  cbind(
    ifelse(ok, tab$estimate - crit * tab$std.error, NA_real_),
    ifelse(ok, tab$estimate + crit * tab$std.error, NA_real_)
  )
}

fc_check_level <- function(level, arg = "level", call = rlang::caller_env()) {
  if (!fc_is_number(level) || level <= 0 || level >= 1) {
    fc_abort("{.arg {arg}} must be a single number between 0 and 1.",
      call = call
    )
  }
}

#' @export
confint.fewclust <- function(object, parm, level = 0.95, ...) {
  fc_check_level(level)
  ci <- fc_ci(object, level)
  a <- (1 - level) / 2
  pct <- format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3)
  dimnames(ci) <- list(object$table$term, paste(pct, "%"))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' Predictions from a fewclust fit
#'
#' Marginal (population-averaged) predictions: the linear predictor on the
#' log-odds scale or the probability of the event.
#'
#' @param object A `"fewclust"` fit.
#' @param newdata Optional data frame with the covariates of the mean model;
#'   by default the observations used in the fit.
#' @param type `"link"` (log odds, default) or `"response"` (probability).
#' @param ... Unused.
#' @return A numeric vector.
#' @examples
#' fit <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
#' nd <- data.frame(time = c(0, 1), arm = 1)
#' predict(fit, newdata = nd, type = "response")
#' @export
predict.fewclust <- function(object, newdata = NULL,
                             type = c("link", "response"), ...) {
  type <- rlang::arg_match(type)
  if (is.null(newdata)) {
    x_mat <- object$model_matrix
    off <- object$offset
  } else {
    if (!is.data.frame(newdata)) {
      fc_abort("{.arg newdata} must be a data frame.")
    }
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, as.data.frame(newdata),
      xlev = object$xlevels,
      na.action = stats::na.pass
    )
    x_mat <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    off <- stats::model.offset(mf)
  }
  eta <- drop(x_mat %*% object$coefficients)
  if (!is.null(off)) eta <- eta + off
  if (type == "response") stats::plogis(eta) else eta
}

#' Tidy, glance and augment a fewclust fit
#'
#' The broom verbs, returning tibbles.
#'
#' @param x A `"fewclust"` fit, or a `"fewclust_flip"` test for `tidy()`.
#' @param conf.int Whether to add the confidence interval.
#' @param conf.level Confidence level; by default that of the fit.
#' @param newdata Optional data frame for `augment()`; by default the
#'   observations used in the fit.
#' @param type.predict Scale of `.fitted` in `augment()`: `"response"`
#'   (probability, default) or `"link"` (log odds).
#' @param ... Unused.
#' @return `tidy()`: one row per coefficient, with `term`, `estimate`,
#'   `std.error`, `statistic`, `df`, `p.value` (and `conf.low`,
#'   `conf.high`).  For a sign-flip test: one row with `term`, `statistic`,
#'   `p.value`, `n_groups`, `n_signs` and `method`.  `glance()`: one row
#'   describing the fit.  `augment()`: the data with `.fitted` (and
#'   `.resid`, the response residual, for the data of the fit).
#' @examples
#' fit <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
#' tidy(fit, conf.int = TRUE)
#' glance(fit)
#' augment(fit)
#' @export
# argument names of broom
# nolint start: object_name_linter.
tidy.fewclust <- function(x, conf.int = FALSE, conf.level = x$conf_level,
                          ...) {
  cols <- c("term", "estimate", "std.error", "statistic", "df", "p.value")
  tab <- x$table[cols]
  if (conf.int) {
    fc_check_level(conf.level, "conf.level")
    ci <- fc_ci(x, conf.level)
    tab$conf.low <- ci[, 1L]
    tab$conf.high <- ci[, 2L]
  }
  tibble::as_tibble(tab)
}
# nolint end

#' @rdname tidy.fewclust
#' @export
glance.fewclust <- function(x, ...) {
  mmm <- x$method == "mmm"
  tibble::tibble(
    method = x$method, corstr = x$corstr, correction = x$correction,
    fallback = !is.null(x$fallback), nobs = x$n_obs, n_groups = x$n_groups,
    n_occasions = x$n_occasions, df = x$df,
    logLik = if (mmm) x$engine$logLik else NA_real_,
    converged = if (mmm) {
      x$engine$code %in% c(1L, 2L)
    } else {
      identical(as.integer(x$engine$geese$error), 0L)
    }
  )
}

#' @rdname tidy.fewclust
#' @export
# argument names of broom
# nolint start: object_name_linter.
augment.fewclust <- function(x, newdata = NULL,
                             type.predict = c("response", "link"), ...) {
  type <- rlang::arg_match(type.predict)
  if (is.null(newdata)) {
    out <- x$data
    out$.fitted <- unname(predict(x, type = type))
    out$.resid <- x$y - unname(predict(x, type = "response"))
  } else {
    out <- tibble::as_tibble(newdata)
    out$.fitted <- unname(predict(x, newdata = newdata, type = type))
  }
  out
}
# nolint end
