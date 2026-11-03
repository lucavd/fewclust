# Single point of parsing: formula + data + group (+ time) -> aligned design.
# Every engine consumes this structure, so outcome validation, NA handling
# and the ordering keys are computed once.

# Resolve a column argument: a bare column name, a string, a one-sided
# formula, or a variable holding one of those.
fc_resolve_column <- function(expr, data, env, what,
                              call = rlang::caller_env()) {
  if (is.null(expr)) {
    return(NULL)
  }
  nm <- if (is.symbol(expr) && as.character(expr) %in% names(data)) {
    as.character(expr)
  } else {
    val <- eval(expr, env)
    if (inherits(val, "formula")) val <- all.vars(val)
    if (!is.character(val) || length(val) != 1L) {
      fc_abort("{.arg {what}} must name exactly one column of {.arg data}.",
        call = call
      )
    }
    val
  }
  if (!nm %in% names(data)) {
    fc_abort("{.arg {what}} column {.val {nm}} not found in {.arg data}.",
      call = call
    )
  }
  nm
}

# Validate a binary outcome; return 0/1 integers plus the label of the event.
fc_binary_outcome <- function(y, call = rlang::caller_env()) {
  if (is.logical(y)) {
    return(list(y = as.integer(y), event = "TRUE"))
  }
  if (is.factor(y) && nlevels(y) == 2L) {
    return(list(y = as.integer(y) - 1L, event = levels(y)[2L]))
  }
  if (is.numeric(y) && all(y %in% c(0, 1))) {
    return(list(y = as.integer(y), event = "1"))
  }
  fc_abort(c(
    "The outcome must be binary.",
    "i" = "Use 0/1 values, a logical, or a factor with 2 levels."
  ), call = call)
}

# Random-effects terms written lme4-style, `(... | g)`, without depending on
# lme4: the bar terms of the right-hand side, and the formula without them.
fc_is_bar <- function(e) {
  is.call(e) && as.character(e[[1L]])[1L] %in% c("|", "||")
}

fc_findbars <- function(formula) {
  out <- list()
  walk <- function(e) {
    if (!is.call(e)) {
      return(invisible())
    }
    if (identical(e[[1L]], as.name("(")) && fc_is_bar(e[[2L]])) {
      out[[length(out) + 1L]] <<- e[[2L]]
      return(invisible())
    }
    for (a in as.list(e)[-1L]) walk(a)
  }
  walk(formula[[length(formula)]])
  out
}

fc_nobars <- function(formula) {
  if (!length(fc_findbars(formula))) {
    return(formula)
  }
  tt <- stats::terms(formula)
  labs <- attr(tt, "term.labels")
  keep <- labs[!grepl("|", labs, fixed = TRUE)]
  if (!length(keep)) keep <- "1"
  if (attr(tt, "intercept") == 0L) keep <- c(keep, "0")
  stats::reformulate(keep,
    response = formula[[2L]],
    env = environment(formula)
  )
}

fc_parse <- function(formula, data, group_expr = NULL, time_expr = NULL,
                     env = parent.frame(), full_rank = TRUE,
                     call = rlang::caller_env()) {
  if (!(inherits(formula, "formula") && length(formula) == 3L)) {
    fc_abort(
      paste0(
        "{.arg formula} must be a two-sided formula, ",
        "such as {.code y ~ time + arm}."
      ),
      call = call
    )
  }
  if (!is.data.frame(data)) {
    fc_abort("{.arg data} must be a data frame.", call = call)
  }
  data <- as.data.frame(data)
  bars <- fc_findbars(formula)
  fixed <- fc_nobars(formula)
  group <- fc_resolve_column(group_expr, data, env, "group", call = call)
  time <- fc_resolve_column(time_expr, data, env, "time", call = call)
  if (length(bars) > 1L) {
    fc_abort("Use a single random-effects term; got {length(bars)}.",
      call = call
    )
  }
  if (length(bars) == 1L) {
    bar_group <- deparse1(bars[[1L]][[3L]])
    if (!is.null(group) && bar_group != group) {
      fc_abort(
        paste0(
          "fewclust covers a single grouping factor; ",
          "got: {group}, {bar_group}."
        ),
        call = call
      )
    }
    group <- bar_group
  }
  if (is.null(group)) {
    fc_abort(
      paste0(
        "There is no grouping variable: ",
        "name the column that identifies the groups."
      ),
      call = call
    )
  }
  slope_vars <- character(0)
  if (length(bars) == 1L) slope_vars <- all.vars(bars[[1L]][[2L]])

  used <- unique(c(all.vars(fixed), slope_vars, group, time))
  missing_cols <- setdiff(used, names(data))
  if (length(missing_cols) > 0L) {
    fc_abort("Variable{?s} not in {.arg data}: {.var {missing_cols}}.",
      call = call
    )
  }
  if (!is.null(time) && !is.numeric(data[[time]])) {
    fc_abort("{.arg time} must be a numeric column.", call = call)
  }
  cc <- stats::complete.cases(data[used])
  dcc <- droplevels(data[cc, used, drop = FALSE])

  q_ranef <- 1L
  if (length(bars) == 1L) {
    ranef <- as.formula(paste("~", deparse1(bars[[1L]][[2L]])))
    q_ranef <- ncol(model.matrix(ranef, dcc))
  }
  mf <- model.frame(fixed, dcc)
  mt <- attr(mf, "terms")
  out <- fc_binary_outcome(model.response(mf), call = call)
  x_mat <- model.matrix(mt, mf)
  qr_x <- qr(x_mat)
  if (full_rank && qr_x$rank < ncol(x_mat)) {
    fc_abort(c(
      paste0(
        "The mean model is rank deficient: ",
        "{.var {colnames(x_mat)[qr_x$pivot[-seq_len(qr_x$rank)]]}} ",
        "{?is/are} collinear with the other columns."
      ),
      "i" = "Drop or recode the collinear terms."
    ), call = call)
  }
  id <- factor(dcc[[group]])
  if (nlevels(id) < 2L) {
    fc_abort("There are fewer than 2 groups in {.var {group}}.", call = call)
  }
  if (!is.null(time) && anyDuplicated(dcc[c(group, time)])) {
    fc_abort(c(
      paste0(
        "{.arg time} repeats within a group: each group must have one ",
        "observation per occasion."
      ),
      "i" = paste0(
        "Omit {.arg time} when the observations of a group are not a ",
        "series."
      )
    ), call = call)
  }

  list(
    formula = formula, fixed = fixed, bars = bars,
    has_bars = length(bars) == 1L, group = group, time = time, data = dcc,
    X = x_mat, y = out$y, event = out$event, id = id,
    slope_vars = slope_vars, q_ranef = q_ranef,
    offset = model.offset(mf), n_dropped = sum(!cc), terms = mt,
    xlevels = stats::.getXlevels(mt, mf),
    contrasts = attr(x_mat, "contrasts")
  )
}
