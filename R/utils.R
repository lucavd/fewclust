# Shared internal helpers.

# Conditions.  Every error, warning and message of the package carries a
# class, so that it can be caught or silenced selectively:
# "fewclust_error", "fewclust_warning", "fewclust_message".  The text is
# interpolated by cli in the caller's environment.
fc_abort <- function(message, ..., call = rlang::caller_env(),
                     .envir = parent.frame()) {
  cli::cli_abort(message, ...,
    class = "fewclust_error", call = call,
    .envir = .envir
  )
}

fc_warn <- function(message, ..., .envir = parent.frame()) {
  cli::cli_warn(message, ..., class = "fewclust_warning", .envir = .envir)
}

fc_inform <- function(message, ..., .envir = parent.frame()) {
  cli::cli_inform(message, ..., class = "fewclust_message", .envir = .envir)
}

# Argument predicates.
fc_is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x)
}

fc_is_whole <- function(x) {
  fc_is_number(x) && x == round(x)
}

fc_is_flag <- function(x) {
  isTRUE(x) || isFALSE(x)
}

# Singular or plural noun for a count.
fc_plural <- function(n, one, many = paste0(one, "s")) {
  if (n == 1L) one else many
}

# A structured failure: the reason travels with the object instead of being
# thrown away, so refusal messages and comparison tables can display it.
fc_failure <- function(reason) {
  structure(list(ok = FALSE, reason = reason), class = "fc_failure")
}

# Evaluate `expr`, turning errors into an fc_failure prefixed by `what` and
# collecting warnings as notes instead of letting them reach the user.
fc_quietly <- function(expr, what) {
  notes <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      fc_failure(paste0(
        what, ": ",
        conditionMessage(e)
      ))
    }),
    warning = function(w) {
      notes <<- c(notes, paste0(what, ": ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  list(value = value, notes = unique(notes))
}

# Standard errors from a variance matrix; a non-positive diagonal gives NA
# instead of a warning.
fc_se <- function(vcov) {
  v <- diag(vcov)
  sqrt(ifelse(is.finite(v) & v > 0, v, NA_real_))
}

# The usability criterion of the study, coefficient by coefficient: finite
# estimate and standard error, |estimate| < 10, 0 < SE < 10 on the logit
# scale.
fc_usable <- function(estimate, se) {
  is.finite(estimate) & is.finite(se) & abs(estimate) < 10 &
    se > 0 & se < 10
}

# Wald-type coefficient table with a t(df) or standard normal reference.
# Coefficients that fail the usability criterion keep their estimate and
# standard error but get no statistic, p-value or interval.
fc_coef_table <- function(term, estimate, se, reference, df, alpha) {
  ok <- fc_usable(estimate, se)
  statistic <- ifelse(ok, estimate / se, NA_real_)
  df <- if (reference == "t") df else Inf
  crit <- qt(1 - alpha / 2, df = df)
  data.frame(
    term = term, estimate = unname(estimate), std.error = unname(se),
    statistic = unname(statistic), df = rep(df, length(estimate)),
    p.value = unname(2 * pt(-abs(statistic), df = df)),
    conf.low = unname(ifelse(ok, estimate - crit * se, NA_real_)),
    conf.high = unname(ifelse(ok, estimate + crit * se, NA_real_)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Engine result for a Wald strategy: a failure when no coefficient is
# usable, otherwise the table with a note naming the unusable terms.
fc_wald_result <- function(term, estimate, se, vcov, model, reference,
                           df, alpha, notes = character(0)) {
  ok <- fc_usable(estimate, se)
  if (!any(ok)) {
    return(fc_failure(paste0(
      "no usable estimate: every coefficient or its standard error is ",
      "non-finite or of 10 or more on the logit scale"
    )))
  }
  if (!all(ok)) {
    notes <- c(notes, sprintf(paste0(
      "no usable inference for %s: estimate or standard error non-finite ",
      "or of 10 or more on the logit scale (statistic, p-value and interval ",
      "set to NA)"
    ), paste(term[!ok], collapse = ", ")))
  }
  list(
    ok = TRUE,
    coefficients = fc_coef_table(
      term, estimate, se, reference, df,
      alpha
    ),
    vcov = vcov, model = model, reference = reference, notes = notes
  )
}
