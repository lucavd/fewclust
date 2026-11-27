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
  tg <- if (is.null(parsed$time)) {
    list(regime = NA_character_, regular = TRUE, groups_with_gaps = 0L)
  } else {
    fc_time_grid(parsed$data[[parsed$time]], id)
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
      time_grid_regime = tg$regime,
      time_grid_regular = tg$regular,
      groups_with_gaps = tg$groups_with_gaps,
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

# Single source of truth for time-grid classification, internal gap counts,
# and the MMM occasion index Z:
#  1. "regular": exact equally spaced shared grid (occasion = 0, 1, ..., K-1,
#     matching MMLB::mm on the study's designs).
#  2. "regular_missing": shared grid whose values are integer multiples of a
#     common step h with >= 75% of slots observed (e.g. 0, 1, 2, 4 missing 3);
#     occasion = (tv - min) / h and gaps are counted on the full regular grid.
#  3. "shared_irregular": shared planned irregular schedule (e.g. 0, 1, 2, 4,
#     8); occasion = (tv - min(grid)) / mean(diff(grid)) (independent of
#     missingness, span K - 1) and gaps are counted on the planned grid.
#  4. "continuous": truly continuous / jittered timestamps; no discrete gaps
#     (0L) and occasion = (tv - min(tv)) / median(within-group diffs).
fc_time_grid <- function(tv, id) {
  grid <- sort(unique(tv))
  n_grid <- length(grid)
  count_gaps <- function(rk) {
    as.integer(sum(tapply(
      rk, id, function(r) length(r) > 1L && any(diff(sort(r)) > 1L)
    )))
  }
  if (n_grid < 3L) {
    return(list(
      regime = "regular", regular = TRUE, shared = TRUE,
      occasion = match(tv, grid) - 1, groups_with_gaps = 0L
    ))
  }
  d <- diff(grid)
  if (isTRUE(all.equal(d, rep(d[1L], length(d))))) {
    rk <- match(tv, grid)
    return(list(
      regime = "regular", regular = TRUE, shared = TRUE,
      occasion = rk - 1, groups_with_gaps = count_gaps(rk)
    ))
  }
  max_n_i <- max(as.integer(table(id)))
  counts <- as.integer(table(tv))
  is_shared <- mean(counts >= 2L) >= 0.5 && n_grid <= 2L * max_n_i
  if (is_shared) {
    h <- min(d)
    mult <- (grid - grid[1L]) / h
    is_mult <- all(abs(mult - round(mult)) < 1e-6)
    n_full <- round(mult[n_grid]) + 1
    if (is_mult && n_grid / n_full >= 0.75) {
      occ <- round((tv - grid[1L]) / h)
      return(list(
        regime = "regular_missing", regular = FALSE, shared = TRUE,
        occasion = occ, groups_with_gaps = count_gaps(occ + 1)
      ))
    }
    rk <- match(tv, grid)
    return(list(
      regime = "shared_irregular", regular = FALSE, shared = TRUE,
      occasion = (tv - grid[1L]) / mean(d), groups_with_gaps = count_gaps(rk)
    ))
  }
  ord <- order(id, tv)
  steps <- unlist(tapply(tv[ord], id[ord], diff), use.names = FALSE)
  step_unit <- stats::median(steps[steps > 0])
  occ <- if (is.finite(step_unit) && step_unit > 0) {
    (tv - min(tv)) / step_unit
  } else {
    stats::ave(seq_along(tv), id, FUN = seq_along) - 1
  }
  list(
    regime = "continuous", regular = FALSE, shared = FALSE,
    occasion = occ, groups_with_gaps = 0L
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
