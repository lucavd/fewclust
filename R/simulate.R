#' Simulate clustered binary outcome data
#'
#' Generates data from a random-intercept (optionally random-slope) logistic
#' model: observations nested in groups, a group-level treatment indicator
#' `arm` (balanced by design), a within-group measurement occasion `time`
#' scaled to \eqn{[0, 1]}, and an observation-level covariate `x`.  The
#' package's examples, tests and vignettes all draw their data from this
#' generator, so everything runs without external data.
#'
#' The intercept is solved from `prevalence` using the logistic-normal
#' attenuation approximation, so the realised marginal prevalence is close to
#' the requested one for moderate `sigma_u`.
#'
#' @param n_groups Number of groups (clusters).  At least 2.
#' @param n_occasions Number of measurement occasions per group.  At least 1.
#' @param effect Log-odds coefficient of the group-level `arm` indicator.
#' @param sigma_u Standard deviation of the random intercept.
#' @param slope_sd Standard deviation of the group-varying slope on `time`;
#'   `0` (the default) generates random intercepts only.
#' @param prevalence Target marginal prevalence of the outcome, in (0, 1).
#' @param seed Optional integer seed.  When supplied, the caller's RNG state
#'   is left untouched.
#'
#' @return A tibble with `n_groups * n_occasions` rows and columns
#'   `id` (factor), `time` (numeric in \eqn{[0, 1]}), `arm` (0/1, constant
#'   within group), `x` (numeric) and `y` (0/1).
#'
#' @examples
#' d <- simulate_clustered_binary(n_groups = 12, n_occasions = 4, seed = 1)
#' d
#' fit_fewclust(d, y ~ time + arm, id = id, time = time)
#' @export
simulate_clustered_binary <- function(n_groups = 24, n_occasions = 4,
                                      effect = 0.6, sigma_u = 1, slope_sd = 0,
                                      prevalence = 0.3, seed = NULL) {
  stopifnot(
    "`n_groups` must be a single number >= 2" =
      is.numeric(n_groups) && length(n_groups) == 1L && n_groups >= 2,
    "`n_occasions` must be a single number >= 1" =
      is.numeric(n_occasions) && length(n_occasions) == 1L && n_occasions >= 1,
    "`effect` must be a single number" =
      is.numeric(effect) && length(effect) == 1L,
    "`sigma_u` must be a single non-negative number" =
      is.numeric(sigma_u) && length(sigma_u) == 1L && sigma_u >= 0,
    "`slope_sd` must be a single non-negative number" =
      is.numeric(slope_sd) && length(slope_sd) == 1L && slope_sd >= 0,
    "`prevalence` must be a single number strictly inside (0, 1)" =
      fc_is_number(prevalence) && prevalence > 0 && prevalence < 1
  )
  n_groups <- as.integer(n_groups)
  n_occasions <- as.integer(n_occasions)

  if (!is.null(seed)) {
    old <- globalenv()$.Random.seed
    on.exit(fc_restore_seed(old), add = TRUE)
    set.seed(seed)
  }

  id <- factor(rep(seq_len(n_groups), each = n_occasions))
  time <- rep(
    if (n_occasions == 1L) 0 else seq(0, 1, length.out = n_occasions),
    times = n_groups
  )
  arm <- rep(rep(c(0L, 1L), length.out = n_groups), each = n_occasions)
  x <- rnorm(n_groups * n_occasions)

  u <- rnorm(n_groups, sd = sigma_u)
  s <- rnorm(n_groups, sd = slope_sd)

  # Logistic-normal attenuation: marginal p ~ plogis(eta / sqrt(1 + 0.346 s^2)).
  intercept <- qlogis(prevalence) * sqrt(1 + 0.346 * sigma_u^2) -
    effect * mean(arm) - 0.4 * mean(time)
  eta <- intercept + effect * arm + 0.4 * time + 0.3 * x +
    u[as.integer(id)] + s[as.integer(id)] * time
  y <- rbinom(length(eta), size = 1L, prob = plogis(eta))

  tibble::tibble(id = id, time = time, arm = arm, x = x, y = y)
}
