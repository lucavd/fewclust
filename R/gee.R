# GEE engine in two layers: fc_gee_model() fits geeglm once and returns the
# model with its notes, or a failure; fc_gee_infer() derives the sandwich
# variance and the coefficient table from it, or a failure.
# The fit is the one of the simulation study: geepack::geeglm with the
# scale parameter fixed at 1 (scale.fix = TRUE) on rows ordered by group
# and, within a group, by `time` when it is given -- with an AR(1) working
# correlation the lags are the consecutive rows of a group.  The response
# is the validated 0/1 vector from fc_parse; warnings from the fit are kept
# as notes.

fc_gee_model <- function(parsed, corstr) {
  ord <- if (is.null(parsed$time)) {
    order(parsed$id)
  } else {
    order(parsed$id, parsed$data[[parsed$time]])
  }
  dat <- parsed$data[ord, , drop = FALSE]
  dat[[".fc_y"]] <- parsed$y[ord]
  dat[[".fc_id"]] <- parsed$id[ord]
  f <- parsed$fixed
  f[[2L]] <- quote(.fc_y)
  q <- fc_quietly(
    geepack::geeglm(
      f,
      id = .fc_id, data = dat, family = binomial("logit"), corstr = corstr,
      scale.fix = TRUE
    ),
    "GEE fit"
  )
  if (inherits(q$value, "fc_failure")) {
    return(q$value)
  }
  err <- q$value$geese$error
  if (!is.null(err) && !is.na(err) && err != 0L) {
    out <- fc_failure(sprintf(
      "GEE fit did not converge (geese error %d)", err
    ))
    out$corstr_degenerate <- corstr != "independence"
    return(out)
  }
  list(model = q$value, notes = q$notes)
}

fc_gee_infer <- function(m, notes, control, type, reference, alpha, df,
                         scale = 1) {
  v <- tryCatch(gee_vcov(m, type, fg_bound = control$fg_bound),
    error = function(e) {
      msg <- conditionMessage(e)
      out <- fc_failure(msg)
      out$corstr_degenerate <- !identical(m$corstr, "independence") &&
        grepl("not positive definite", msg, fixed = TRUE)
      out
    }
  )
  if (inherits(v, "fc_failure")) {
    return(v)
  }
  fc_wald_result(
    names(coef(m)), coef(m), fc_se(v), v, m, reference, df,
    alpha, notes,
    scale = scale
  )
}

fit_gee <- function(parsed, control, type, reference, alpha, df,
                    corstr = "ar1") {
  base <- fc_gee_model(parsed, corstr)
  if (inherits(base, "fc_failure")) {
    return(base)
  }
  fc_gee_infer(
    base$model, base$notes, control, type, reference, alpha, df,
    scale = fc_col_scale(parsed$X)
  )
}
