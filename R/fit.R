#' Fit a model for a clustered binary outcome
#'
#' Fits a regression model for a binary outcome measured repeatedly within a
#' small number of groups (clusters), with tests and intervals that remain
#' calibrated with few groups.  The default is the method recommended by
#' the simulation study behind the package: the marginalized multilevel
#' model of Heagerty (1999) with a convergence check, its Wald statistics
#' referred to a t distribution on N - p degrees of freedom (N groups, p
#' coefficients).  When that model gives no usable estimate, GEE with an
#' AR(1) working correlation and the Fay-Graubard small-sample correction,
#' with the same t reference, is fitted instead, and a message says so.
#'
#' The coefficients are marginal (population-averaged) log odds ratios.
#'
#' @section Methods:
#' * `"mmm"`: the marginalized multilevel model (Heagerty 1999).  With
#'   `time` the latent standard deviation follows a linear trend over the
#'   occasions and a first-order transition term links consecutive
#'   observations; the convergence check (see [fewclust_control()])
#'   replaces a trend fit at the boundary with the fit with a constant
#'   latent standard deviation when that fit converges and has finite
#'   estimates and standard errors; otherwise the model gives no result
#'   (`"auto"` then fits GEE, `"mmm"` stops with an error).  The fit at the
#'   boundary is never reported.  Without `time` the latent standard
#'   deviation is constant and there is no transition term.  Standard
#'   errors are model-based.
#' * `"gee"`: [geepack::geeglm()] with the scale parameter fixed at 1, on
#'   rows ordered by group and `time`, with the small-sample sandwich
#'   variance chosen by `correction` (see [gee_vcov()]).
#' * `"auto"` (default): `"mmm"`, and `"gee"` (with `corstr` and
#'   `correction`) when `"mmm"` gives no usable estimate.
#'
#' A coefficient is usable when its estimate and standard error are finite
#' and smaller than 10 in absolute value on the logit scale; otherwise its
#' statistic, p-value and interval are `NA` (with a warning).  A model with
#' no usable coefficient gives no result.
#'
#' @section Errors, warnings, messages and notes:
#' The function stops with an error (class `"fewclust_error"`) when no
#' estimate can be computed: an outcome that is not binary or is constant,
#' an outcome that varies within one group only, separation, as many
#' coefficients as groups or more, a single observation per group, repeated
#' values of `time` within a group, random-effects terms in the formula,
#' or no usable estimate from the requested method(s).  It warns (class
#' `"fewclust_warning"`) when the result is less reliable: fewer than 10
#' groups; GEE with 2 occasions and fewer than 15 groups; a GEE correction
#' other than Fay-Graubard; coefficients without inference; a marginalized
#' model that did not converge; series with gaps or unequally spaced
#' occasions.  Messages (class `"fewclust_message"`) report a fallback to
#' GEE and a missing `time`.  Notes stored in the fit and shown by
#' `print()` describe small designs and the convergence check.
#'
#' @param data A data frame or tibble with one row per observation.
#' @param formula A two-sided formula for the mean model, for example
#'   `y ~ time * arm + x`.  The outcome must be binary: 0/1, logical, or a
#'   factor with two levels (the second level is the event).
#' @param id The column that identifies the groups (clusters).  A bare
#'   column name or a string; inside your own functions use `{{ }}`.
#' @param time Optional numeric column with the measurement occasion.  It
#'   orders the observations within a group and gives the marginalized
#'   model its serial dependence; recommended for repeated measurements.
#' @param method `"auto"` (default), `"mmm"` or `"gee"`.
#' @param corstr Working correlation of GEE: `"ar1"` (default),
#'   `"exchangeable"` or `"independence"`.
#' @param correction Small-sample correction of the GEE sandwich variance:
#'   `"fg"` (Fay-Graubard, default), `"md"` (Mancl-DeRouen), `"kc"`
#'   (Kauermann-Carroll) or `"none"`.
#' @param conf_level Confidence level of the intervals.
#' @param control Tuning constants from [fewclust_control()].
#'
#' @return An object of class `"fewclust"`.  Use `print()`, `summary()`,
#'   `coef()`, `vcov()`, `confint()`, `predict()`, `nobs()`, and the broom
#'   verbs [tidy()], [glance()] and [augment()], which return tibbles.  The
#'   list also holds `method` (`"mmm"` or `"gee"`, the model reported),
#'   `fallback` (`NULL`, or why the marginalized model was not used),
#'   `df`, `notes`, `convergence` (the record of the convergence check of
#'   the marginalized model, also when it gave no result and GEE was
#'   reported instead) and
#'   `engine` (the underlying fit).
#'
#' @references
#' Heagerty PJ (1999). Marginally specified logistic-normal models for
#' longitudinal binary data. *Biometrics* 55(3), 688-698.
#' \doi{10.1111/j.0006-341X.1999.00688.x}
#'
#' Fay MP, Graubard BI (2001). Small-sample adjustments for Wald-type tests
#' using sandwich estimators. *Biometrics* 57(4), 1198-1206.
#' \doi{10.1111/j.0006-341X.2001.01198.x}
#'
#' @examples
#' fit <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
#' fit
#' tidy(fit, conf.int = TRUE)
#'
#' # GEE with an exchangeable working correlation
#' fit_fewclust(fc_trial, y ~ time * arm,
#'   id = id, time = time,
#'   method = "gee", corstr = "exchangeable"
#' )
#'
#' # in a pipe
#' fc_trial |>
#'   fit_fewclust(y ~ time + arm + x, id = id, time = time) |>
#'   glance()
#' @seealso [flip_test()] for a p-value with very small designs,
#'   [gee_vcov()] for the sandwich estimators.
#' @export
fit_fewclust <- function(data, formula, id, time = NULL,
                         method = c("auto", "mmm", "gee"),
                         corstr = c("ar1", "exchangeable", "independence"),
                         correction = c("fg", "md", "kc", "none"),
                         conf_level = 0.95, control = fewclust_control()) {
  cl <- match.call()
  call <- rlang::current_env()
  method <- rlang::arg_match(method)
  corstr <- rlang::arg_match(corstr)
  correction <- rlang::arg_match(correction)
  if (!fc_is_number(conf_level) || conf_level <= 0 || conf_level >= 1) {
    fc_abort(paste0(
      "{.arg conf_level} must be a single number between 0 and 1, ",
      "such as 0.95."
    ))
  }
  if (!inherits(control, "fewclust_control")) {
    fc_abort("{.arg control} must be created with {.fn fewclust_control}.")
  }
  cols <- fc_columns(data, formula, rlang::enquo(id), rlang::enquo(time),
    id_missing = missing(id), call = call
  )
  if (length(fc_findbars(formula))) {
    fc_abort(c(
      paste0(
        "{.arg formula} must not contain random-effects terms ",
        "such as {.code (1 | id)}."
      ),
      "i" = paste0(
        "Give the groups with {.arg id}: the methods model the ",
        "dependence within a group."
      )
    ))
  }
  parsed <- fc_parse(formula, data, cols$id, cols$time, call = call)
  dg <- fc_diagnose(parsed)
  fc_check_design(dg, call = call)
  if (is.null(parsed$time)) fc_inform_no_time()

  alpha <- 1 - conf_level
  type <- c(fg = "fg", md = "md", kc = "kc", none = "unc")[[correction]]
  fit_gee_t <- function() {
    fit_gee(parsed, control, type, "t", alpha, dg$df_t, corstr = corstr)
  }
  fallback <- NULL
  mmm_check <- NULL
  if (method == "gee") {
    res <- fit_gee_t()
    if (inherits(res, "fc_failure")) {
      fc_abort(c(
        "GEE gave no usable estimate: {res$reason}.",
        "i" = paste0(
          "Try the default {.code method = \"auto\"}, or {.fn flip_test} ",
          "for a p-value."
        )
      ))
    }
  } else {
    res <- fit_mmm(parsed, control, "t", alpha, dg$df_t)
    if (inherits(res, "fc_failure") && method == "mmm") {
      fc_abort(c(
        paste0(
          "The marginalized multilevel model gave no usable estimate: ",
          "{res$reason}."
        ),
        "i" = paste0(
          "Try {.code method = \"gee\"}, or {.fn flip_test} for a p-value."
        )
      ), convergence = res$check)
    }
    if (inherits(res, "fc_failure")) {
      fallback <- res$reason
      mmm_check <- res$check
      res <- fit_gee_t()
      if (inherits(res, "fc_failure")) {
        fc_abort(c(
          "Neither model gave a usable estimate.",
          "x" = "Marginalized multilevel model: {fallback}.",
          "x" = "GEE: {res$reason}.",
          "i" = paste0(
            "{.fn flip_test} gives a p-value for one coefficient without ",
            "fitting either model."
          )
        ), convergence = mmm_check)
      }
      fc_inform(c(
        paste0(
          "The marginalized multilevel model gave no usable estimate ",
          "({fallback})."
        ),
        "i" = paste0(
          "Reporting GEE with {fc_gee_label(corstr, correction)} instead."
        )
      ))
    }
  }
  out <- fc_new_fit(res, parsed, dg, fallback, corstr, correction, conf_level,
    control,
    call = cl, mmm_check = mmm_check
  )
  fc_fit_warnings(out)
  out$notes <- fc_fit_notes(out, res$notes)
  out
}

fc_new_fit <- function(res, parsed, dg, fallback, corstr, correction,
                       conf_level, control, call, mmm_check = NULL) {
  used <- if (inherits(res$model, "geeglm")) "gee" else "mmm"
  tab <- tibble::as_tibble(res$coefficients)
  vc <- res$vcov
  dimnames(vc) <- list(tab$term, tab$term)
  structure(
    list(
      coefficients = stats::setNames(tab$estimate, tab$term),
      vcov = vc, table = tab, method = used, fallback = fallback,
      corstr = if (used == "gee") corstr else NA_character_,
      correction = if (used == "gee") correction else NA_character_,
      df = dg$df_t, conf_level = conf_level, engine = res$model,
      convergence = if (used == "mmm") res$model$check else mmm_check,
      n_obs = dg$n_obs, n_groups = dg$n_groups,
      n_occasions = dg$n_occasions, notes = character(0), call = call,
      formula = parsed$formula, terms = parsed$terms,
      xlevels = parsed$xlevels, contrasts = parsed$contrasts,
      model_matrix = parsed$X, offset = parsed$offset, y = parsed$y,
      data = tibble::as_tibble(parsed$data), id = parsed$group,
      time = parsed$time, event = parsed$event, design = dg,
      control = control
    ),
    class = "fewclust"
  )
}

fc_inform_no_time <- function() {
  fc_inform(
    c(
      paste0(
        "No {.arg time} given: the rows of each group are taken in the ",
        "order of {.arg data}."
      ),
      "i" = paste0(
        "For repeated measurements give {.arg time}, so that the series ",
        "are ordered and their serial dependence is modelled."
      )
    ),
    .frequency = "once", .frequency_id = "fewclust_no_time"
  )
}

# Resolve `data`, `formula`, `id` and `time` (tidyselect) to column names.
fc_columns <- function(data, formula, id, time, id_missing,
                       call = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    fc_abort("{.arg data} must be a data frame or a tibble.", call = call)
  }
  if (!(inherits(formula, "formula") && length(formula) == 3L)) {
    fc_abort(
      paste0(
        "{.arg formula} must be a two-sided formula, ",
        "such as {.code y ~ time + arm}."
      ),
      call = call
    )
  }
  if (id_missing) {
    fc_abort(c(
      "{.arg id} is missing.",
      "i" = paste0(
        "Name the column that identifies the groups, ",
        "for example {.code id = patient}."
      )
    ), call = call)
  }
  time_col <- NULL
  if (!rlang::quo_is_null(time)) {
    time_col <- fc_select_column(time, data, "time", call)
  }
  list(id = fc_select_column(id, data, "id", call), time = time_col)
}

fc_select_column <- function(quo, data, what, call) {
  pos <- tidyselect::eval_select(quo, data,
    allow_rename = FALSE,
    error_call = call
  )
  if (length(pos) != 1L) {
    fc_abort("{.arg {what}} must select one column, not {length(pos)}.",
      call = call
    )
  }
  names(pos)
}

fc_gee_label <- function(corstr, correction) {
  cor <- c(
    ar1 = "an AR(1)", exchangeable = "an exchangeable",
    independence = "an independence"
  )
  corr <- c(
    fg = "the Fay-Graubard correction", md = "the Mancl-DeRouen correction",
    kc = "the Kauermann-Carroll correction",
    none = "no small-sample correction"
  )
  paste(cor[[corstr]], "working correlation and", corr[[correction]])
}

# Errors: designs in which no estimate can be computed.
fc_check_design <- function(dg, call = rlang::caller_env()) {
  if (dg$df_t < 1L) {
    fc_abort(c(
      paste0(
        "The mean model has {dg$p_fixed} coefficients for {dg$n_groups} ",
        "groups: there must be more groups than coefficients."
      ),
      "i" = "Use a smaller mean model."
    ), call = call)
  }
  if (dg$group_sizes[["max"]] == 1) {
    fc_abort(c(
      paste0(
        "Every group has a single observation, so there is no clustering ",
        "to account for."
      ),
      "i" = paste0(
        "Use an ordinary logistic regression, {.fn glm} with ",
        "{.code family = binomial}."
      )
    ), call = call)
  }
  if (dg$prevalence %in% c(0, 1)) {
    fc_abort(c(
      paste0(
        "The outcome is constant (always {dg$prevalence}): ",
        "there is nothing to model."
      ),
      "i" = "Check the outcome column and the rows with missing values."
    ), call = call)
  }
  mixed <- dg$n_groups - dg$pure_groups
  if (mixed < 2L) {
    fc_abort(c(
      paste0(
        if (mixed == 1L) {
          "The outcome varies in only one group: "
        } else {
          "The outcome does not vary within any group: "
        },
        "every other group has a constant outcome."
      ),
      "i" = paste0(
        "These models need the outcome to vary within at least two ",
        "groups; consider exact or permutation inference at the group level."
      )
    ), call = call)
  }
  if (dg$separation$status != "none") {
    ev <- gsub("{", "{{", dg$separation$evidence, fixed = TRUE)
    ev <- gsub("}", "}}", ev, fixed = TRUE)
    fc_abort(c(
      paste0(
        "The mean model shows {dg$separation$status} separation: some ",
        "covariate predicts the outcome perfectly, so the estimates are ",
        "infinite."
      ),
      stats::setNames(ev, rep("*", length(ev))),
      "i" = paste0(
        "Remove or coarsen the covariate, or use a penalised (Firth) ",
        "logistic regression, as in the packages {.pkg logistf} or ",
        "{.pkg brglm2}."
      )
    ), call = call)
  }
  invisible(dg)
}

# Warnings: results that are less reliable.
fc_fit_warnings <- function(fit) {
  dg <- fit$design
  gee <- fit$method == "gee"
  if (dg$n_groups < 10L) {
    fc_warn(c(
      paste0(
        "Only {dg$n_groups} groups: with fewer than 10 the tests and ",
        "intervals may be unreliable."
      ),
      "i" = paste0(
        "{.fn flip_test} gives a p-value that keeps its level with few ",
        "groups."
      )
    ))
  }
  if (gee && dg$n_occasions <= 2 && dg$n_groups < 15L) {
    fc_warn(c(
      paste0(
        "GEE with 2 occasions and fewer than 15 groups tends to reject ",
        "too often."
      ),
      "i" = "Check the result with {.fn flip_test}."
    ))
  }
  if (gee && fit$correction == "none") {
    fc_warn(c(
      paste0(
        "GEE without a small-sample correction: its standard errors are ",
        "too small with few groups."
      ),
      "i" = "Use the default {.code correction = \"fg\"} (Fay-Graubard)."
    ))
  }
  if (gee && fit$correction %in% c("md", "kc")) {
    fc_warn(c(
      paste0(
        "GEE with the ",
        c(md = "Mancl-DeRouen", kc = "Kauermann-Carroll")[[fit$correction]],
        " correction."
      ),
      "i" = paste0(
        "With few groups the recommended correction is the default, ",
        "{.code correction = \"fg\"} (Fay-Graubard)."
      )
    ))
  }
  fc_series_warnings(fit, dg, gee)
  bad <- fit$table$term[is.na(fit$table$statistic)]
  if (length(bad)) {
    fc_warn(c(
      paste0(
        "No inference for {.var {bad}}: the estimate or its standard error ",
        "is not finite or is 10 or more on the logit scale."
      ),
      "i" = paste0(
        "The statistic, p-value and interval are {.val {NA}}; check the ",
        "scale of the covariate and for sparse cells."
      )
    ))
  }
  if (!gee && !fit$engine$code %in% c(1L, 2L)) {
    fc_warn(c(
      paste0(
        "The optimiser of the marginalized multilevel model did not ",
        "converge (nlm code {fit$engine$code})."
      ),
      "i" = paste0(
        "Treat the estimates with caution, or compare with ",
        "{.code method = \"gee\"}."
      )
    ))
  }
  invisible(fit)
}

fc_series_warnings <- function(fit, dg, gee) {
  if (is.null(fit$time) || (gee && fit$corstr != "ar1")) {
    return(invisible(fit))
  }
  if (isTRUE(dg$groups_with_gaps > 0L)) {
    fc_warn(c(
      paste0(
        "Some series miss an occasion ({dg$groups_with_gaps} of ",
        "{dg$n_groups} groups)."
      ),
      "i" = paste0(
        "Consecutive observations of a group are treated as one occasion ",
        "apart."
      )
    ))
  }
  if (!dg$time_grid_regular) {
    fc_warn(c(
      "The values of {.arg time} are not equally spaced.",
      "i" = if (gee) {
        paste0(
          "The AR(1) working correlation treats consecutive observations ",
          "as one lag apart."
        )
      } else {
        paste0(
          "The marginalized multilevel model uses the rank of each value ",
          "as the occasion index."
        )
      }
    ))
  }
  invisible(fit)
}

# Notes shown by print(): facts about the fit that need no action.
fc_fit_notes <- function(fit, engine_notes) {
  dg <- fit$design
  as_warning <- startsWith(engine_notes, "no usable inference for ")
  c(
    if (dg$n_dropped > 0L) {
      sprintf(
        "%d %s with missing values %s dropped.", dg$n_dropped,
        fc_plural(dg$n_dropped, "row"), fc_plural(dg$n_dropped, "was", "were")
      )
    },
    if (dg$n_occasions <= 3 && dg$n_groups <= 25L) {
      sprintf(
        paste0(
          "With %g occasions and %d groups the tests have little power and ",
          "some datasets give no usable estimate; flip_test() gives a ",
          "p-value that keeps its level."
        ),
        dg$n_occasions, dg$n_groups
      )
    },
    if (fit$method == "mmm" && dg$n_groups <= 10L) {
      paste0(
        "With 10 groups or fewer the marginalized multilevel model tends ",
        "to be conservative: its tests reject less often than the nominal ",
        "level."
      )
    },
    unique(engine_notes[!as_warning])
  )
}
