# Measurements of the design (internal).  They drive the errors, warnings
# and notes of fit_fewclust(): number of groups and occasions, prevalence,
# groups with a constant outcome, separation, gaps in the series.

fc_diagnose <- function(parsed) {
  y <- parsed$y
  id <- parsed$id
  sizes <- as.numeric(table(id))
  n_groups <- nlevels(id)
  distinct_per_group <- function(v, summary) {
    summary(tapply(v, id, function(z) length(unique(z))))
  }
  # T of the study is the design maximum: the longest series observed
  n_occasions <- if (is.null(parsed$time)) {
    max(sizes)
  } else {
    distinct_per_group(parsed$data[[parsed$time]], max)
  }
  occ_slope <- if (length(parsed$slope_vars)) {
    min(vapply(parsed$slope_vars, function(v) {
      distinct_per_group(parsed$data[[v]], stats::median)
    }, numeric(1)))
  } else {
    NA_real_
  }
  prev_by_group <- tapply(y, id, mean)
  grid <- NA
  if (!is.null(parsed$time)) grid <- sort(unique(parsed$data[[parsed$time]]))
  # a gap is a missing occasion of the grid between the first and the last
  # occasion of a group (a series truncated at either end has none)
  gaps <- if (is.null(parsed$time)) {
    0L
  } else {
    rk <- match(parsed$data[[parsed$time]], grid)
    sum(tapply(rk, id, function(r) length(r) > 1L && any(diff(sort(r)) > 1L)))
  }
  structure(
    list(
      n_obs = length(y), n_groups = n_groups, group = parsed$group,
      time = parsed$time, event = parsed$event,
      group_sizes = c(
        min = min(sizes), median = stats::median(sizes),
        max = max(sizes)
      ),
      n_occasions = n_occasions,
      # the occasion index of the dependence model is the rank of `time` in
      # this grid; it matches the study's occasions only when the grid is
      # equally spaced (one grid shared by all groups)
      time_grid_regular = is.null(parsed$time) || length(grid) < 3L ||
        isTRUE(all.equal(diff(grid), rep(diff(grid)[1L], length(grid) - 1L))),
      groups_with_gaps = as.integer(gaps),
      prevalence = mean(y),
      pure_groups = sum(prev_by_group %in% c(0, 1)),
      icc = fc_icc(y, id),
      separation = detect_separation(parsed$X, y),
      slope_vars = parsed$slope_vars, q_ranef = parsed$q_ranef,
      # the random-effects covariance needs more distinct occasions of the
      # slope variables than random-effect terms (two subject-level
      # parameters saturate two observations)
      ranef_identifiable = is.na(occ_slope) || occ_slope > parsed$q_ranef,
      p_fixed = ncol(parsed$X),
      df_t = n_groups - ncol(parsed$X),
      n_dropped = parsed$n_dropped
    ),
    class = "cluster_diagnosis"
  )
}

# ANOVA (moment) estimator of the intraclass correlation on the 0/1 scale.
fc_icc <- function(y, id) {
  n <- length(y)
  g <- nlevels(id)
  ni <- as.numeric(table(id))
  if (n <= g) {
    return(NA_real_)
  }
  ybar_i <- tapply(y, id, mean)
  msb <- sum(ni * (ybar_i - mean(y))^2) / (g - 1)
  msw <- sum((y - ybar_i[as.integer(id)])^2) / (n - g)
  n0 <- (n - sum(ni^2) / n) / (g - 1)
  denom <- msb + (n0 - 1) * msw
  if (denom == 0) 0 else (msb - msw) / denom
}


#' @export
print.cluster_diagnosis <- function(x, ...) {
  cat("<cluster_diagnosis>\n")
  cat(sprintf(
    "  observations : %d in %d groups (sizes %g-%g, median %g)\n",
    x$n_obs, x$n_groups, x$group_sizes[["min"]],
    x$group_sizes[["max"]], x$group_sizes[["median"]]
  ))
  cat(sprintf(
    "  occasions    : %g at most (%s)\n", x$n_occasions,
    if (is.null(x$time)) {
      "largest group; no `time` given"
    } else {
      paste0("distinct values of ", x$time)
    }
  ))
  cat(sprintf(
    "  outcome      : event \"%s\", prevalence %.3f, %d pure group%s\n",
    x$event, x$prevalence, x$pure_groups,
    if (x$pure_groups == 1L) "" else "s"
  ))
  cat(sprintf("  dependence   : ANOVA ICC %.3f\n", x$icc))
  cat(sprintf(
    "  separation   : %s%s\n", x$separation$status,
    if (length(x$separation$evidence)) {
      paste0(" (", x$separation$evidence[1L], ")")
    } else {
      ""
    }
  ))
  cat(sprintf(
    "  mean model   : %d coefficients; t reference on %d df\n",
    x$p_fixed, x$df_t
  ))
  if (length(x$slope_vars)) {
    cat(sprintf(
      "  random part  : %d terms (slopes: %s); covariance identifiable: %s\n",
      x$q_ranef, paste(x$slope_vars, collapse = ", "),
      if (x$ranef_identifiable) "yes" else "NO"
    ))
  }
  if (x$n_dropped > 0L) {
    cat(sprintf(
      "  note         : %d incomplete row%s dropped\n", x$n_dropped,
      if (x$n_dropped == 1L) "" else "s"
    ))
  }
  invisible(x)
}
