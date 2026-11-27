mmm_design <- function(d, formula, trend, transition) {
  d <- d[order(d$id, if (transition) d$time else seq_len(nrow(d))), ]
  one <- matrix(1, nrow(d), 1L, dimnames = list(NULL, "(Intercept)"))
  fewclust:::fc_mmm_design(
    d$y, model.matrix(formula, d), factor(d$id),
    Z = if (trend) cbind(one, time = d$time) else one,
    W = if (transition) one else NULL, q = 20L)
}

test_that("Gauss-Hermite rule integrates moments of the standard normal", {
  gh <- fewclust:::fc_gh(20L)
  expect_equal(sum(gh$w), 1)
  expect_equal(sum(gh$w * gh$z^2), 1)
  expect_equal(sum(gh$w * gh$z^4), 3)
})

test_that("the marginalized model reaches the optimum of MMLB::mm along the same path", {
  fx <- readRDS(test_path("fixtures", "mmlb_reference.rds"))
  for (nm in names(fx$specs)) {
    sp <- fx$specs[[nm]]
    des <- mmm_design(fx$data[[sp$data]], sp$formula, trend = sp$trend,
                      transition = nm != "clusters")
    fit <- fewclust:::fc_mmm_fit_once(des, fewclust:::fc_mmm_inits(des))
    ref <- fx$reference[fx$reference$fixture == nm, ]
    # columns on the scale of the study are not rescaled, so nlm follows the
    # path of MMLB::mm: observed agreement 1.6e-10 (estimates, relative),
    # 2.4e-7 (standard errors), 7.4e-12 (log-likelihood); tolerances are
    # at least 40 times the observed differences
    expect_equal(unname(fit$theta), ref$estimate, tolerance = 1e-8, info = nm)
    expect_equal(unname(sqrt(diag(fit$vcov))), ref$se, tolerance = 1e-5, info = nm)
    expect_equal(fit$logLik, ref$logLik[1L], tolerance = 1e-9, info = nm)
    expect_equal(fit$code, ref$code[1L], info = nm)
  }
})

test_that("the analytic gradient matches central differences", {
  fx <- readRDS(test_path("fixtures", "mmlb_reference.rds"))
  des <- mmm_design(fx$data$covariate, y ~ time * trt + x, TRUE, TRUE)
  th <- fewclust:::fc_mmm_inits(des) + 0.1
  f <- fewclust:::fc_mmm_negloglik(th, des)
  num <- vapply(seq_along(th), function(j) {
    e <- replace(numeric(length(th)), j, 1e-6)
    (as.numeric(fewclust:::fc_mmm_negloglik(th + e, des)) -
       as.numeric(fewclust:::fc_mmm_negloglik(th - e, des))) / 2e-6
  }, numeric(1))
  expect_equal(unname(attr(f, "gradient")), num, tolerance = 1e-6)
})

test_that("with a vanishing latent variance the likelihood is the logistic GLM", {
  d <- fc_test_data(20, seed = 3)
  des <- mmm_design(d, y ~ time + arm, trend = FALSE, transition = FALSE)
  g <- glm(y ~ time + arm, binomial, d)
  ll <- -as.numeric(fewclust:::fc_mmm_negloglik(c(coef(g) * des$scale, -15), des))
  expect_equal(ll, as.numeric(logLik(g)), tolerance = 1e-8)
})

test_that("the convergence check replaces a boundary fit by the constant-SD fit", {
  d <- fc_test_data(24, seed = 2026)   # the fc_trial generator: no latent trend
  parsed <- fewclust:::fc_parse(y ~ time * arm + x, d, quote(id), quote(time))
  res <- fewclust:::fit_mmm(parsed, fewclust_control(), "t", 0.05, 19L)
  expect_true(res$ok)
  expect_true(res$model$check$requested)
  expect_true(res$model$check$boundary)
  expect_equal(res$model$check$latent_sd, "constant")
  expect_match(res$notes, "convergence check")
  expect_false("log(sigma):occasion" %in% names(res$model$theta))
  # without the check the boundary fit is kept and the switch-off is noted
  raw <- fewclust:::fit_mmm(parsed, fewclust_control(check = FALSE), "t", 0.05,
                            19L)
  expect_true("log(sigma):occasion" %in% names(raw$model$theta))
  expect_match(paste(raw$notes, collapse = " "), "switched off")
  # a steep genuine trend passes the check
  fx <- readRDS(test_path("fixtures", "mmlb_reference.rds"))
  ps <- fewclust:::fc_parse(y ~ time * trt, fx$data$steep, quote(id), quote(time))
  st <- fewclust:::fit_mmm(ps, fewclust_control(), "t", 0.05, 36L)
  expect_false(st$model$check$boundary)
  expect_equal(st$model$check$latent_sd, "time trend")
})

test_that("without `time` the model has a constant latent SD and no transition", {
  d <- fc_test_data(20, seed = 4)
  parsed <- fewclust:::fc_parse(y ~ time + arm, d, quote(id))
  res <- fewclust:::fit_mmm(parsed, fewclust_control(), "t", 0.05, 17L)
  expect_true(res$ok)
  expect_equal(names(res$model$theta),
               c("(Intercept)", "time", "arm", "log(sigma):(Intercept)"))
  expect_false(res$model$check$requested)
  expect_equal(res$coefficients$df, rep(17, 3))
})

test_that("the unit of time changes neither the dependence model nor the optimum", {
  d <- fc_test_data(30, seed = 10)
  days <- d
  days$time <- d$time * 365
  a <- fit_fewclust(d, y ~ time + arm, id = id, time = time, method = "mmm")
  b <- fit_fewclust(days, y ~ time + arm, id = id, time = time, method = "mmm")
  expect_equal(b$method, "mmm")
  keys <- c("requested", "boundary", "latent_sd", "path")
  expect_equal(a$convergence[keys], b$convergence[keys])
  expect_equal(a$engine$logLik, b$engine$logLik, tolerance = 1e-7)
  # the days column is rescaled by 100 for the optimiser, so the two fits
  # reach the optimum along different paths: agreement at nlm's tolerance
  expect_equal(unname(coef(a)[-2]), unname(coef(b)[-2]), tolerance = 1e-4)
  expect_equal(unname(coef(a)[2]), 365 * unname(coef(b)[2]), tolerance = 1e-4)
  expect_equal(tidy(a)$p.value, tidy(b)$p.value, tolerance = 1e-4)
})

test_that("the convergence check catches a false optimum with an indefinite Hessian", {
  d <- simulate_clustered_binary(n_groups = 10, n_occasions = 6, seed = 8)
  d$time <- round(d$time * 5)
  fit <- fit_fewclust(d, y ~ time * arm, id = id, time = time, method = "mmm")
  expect_equal(fit$method, "mmm")
  expect_true(fit$convergence$boundary)
  expect_equal(fit$convergence$latent_sd, "constant")
  expect_match(paste(fit$convergence$reason, collapse = " "), "false optimum")
  # Even with check_coherence = FALSE, the indefinite Hessian at the false
  # optimum yields non-finite SEs and triggers the fallback to constant SD
  raw <- fit_fewclust(d, y ~ time * arm, id = id, time = time, method = "mmm",
                      control = fewclust_control(check_coherence = FALSE))
  expect_true(raw$convergence$boundary)
  expect_equal(raw$convergence$latent_sd, "constant")
  expect_match(paste(raw$convergence$reason, collapse = " "), "not finite")
})

test_that("a trend fit below its own sub-model is treated as a false optimum (stub)", {
  fx <- readRDS(test_path("fixtures", "mmlb_reference.rds"))
  parsed <- fewclust:::fc_parse(y ~ time * trt, fx$data$steep, quote(id), quote(time))
  real <- fewclust:::fc_mmm_fit_once
  local_mocked_bindings(
    fc_mmm_fit_once = function(d, inits) {
      f <- real(d, inits)
      if (d$m_s > 1L) f$logLik <- f$logLik - 5   # pretend the trend fit stopped short
      f
    },
    .package = "fewclust")
  res <- fewclust:::fit_mmm(parsed, fewclust_control(), "t", 0.05, 36L)
  expect_true(res$model$check$boundary)
  expect_equal(res$model$check$latent_sd, "constant")
  expect_match(paste(res$notes, collapse = " "), "false optimum")
  # With check_coherence = FALSE the positive-definite trend fit is kept
  raw <- fewclust:::fit_mmm(
    parsed, fewclust_control(check_coherence = FALSE), "t", 0.05, 36L
  )
  expect_false(raw$model$check$boundary)
  expect_equal(raw$model$check$latent_sd, "time trend")
})

test_that("repeated time values within a group are refused at parse time", {
  d <- fc_test_data(20, seed = 1)
  d$time[d$time == 1] <- 2 / 3
  expect_error(fit_fewclust(d, y ~ time + arm, id = id, time = time),
               "repeats within a group", class = "fewclust_error")
})

test_that("a fit with no usable coefficient is a structured failure", {
  fx <- readRDS(test_path("fixtures", "mmlb_reference.rds"))
  parsed <- fewclust:::fc_parse(y ~ trt + x, fx$data$clusters, quote(id))
  local_mocked_bindings(
    fc_mmm_fit_once = function(d, inits) list(
      theta = setNames(rep(20, 4), c("(Intercept)", "trt", "x", "log(sigma):(Intercept)")),
      vcov = diag(4), logLik = 0, code = 1L, iterations = 1L, p = 3L,
      slope_min = NA_real_),
    .package = "fewclust")
  res <- fewclust:::fit_mmm(parsed, fewclust_control(), "t", 0.05, 37L)
  expect_s3_class(res, "fc_failure")
  expect_match(res$reason, "no usable estimate")
})

test_that("continuous/jittered visit times scale occasion cleanly and fit genuine slopes", {
  set.seed(42)
  d <- simulate_clustered_binary(
    n_groups = 40, n_occasions = 6, slope_sd = 2, sigma_u = 1, seed = 1
  )
  d$day <- d$time * 100 + runif(nrow(d), -2, 2)
  d <- d[order(d$id, d$day), ]
  f_reg <- fit_fewclust(d, y ~ time + arm, id = id, time = time, method = "mmm")
  f_jit <- suppressWarnings(
    fit_fewclust(d, y ~ day + arm, id = id, time = day, method = "mmm")
  )
  expect_equal(f_jit$convergence$path, "trend fit kept")
  expect_true(f_jit$convergence$trend$code %in% c(1L, 2L))
  expect_gt(f_jit$convergence$trend$slope_min, 0.1)
  expect_equal(f_jit$engine$logLik, f_reg$engine$logLik, tolerance = 2)
  # Irregular continuous times should not trigger bogus "miss an occasion" gaps
  expect_equal(f_jit$design$groups_with_gaps, 0L)
})

test_that("many unique jittered times do not cause false convergence or underflow", {
  set.seed(1)
  d <- simulate_clustered_binary(n_groups = 120, n_occasions = 10, seed = 1)
  d$day <- d$time * 100 + runif(nrow(d), 0, 0.5)
  d <- d[order(d$id, d$day), ]
  # Must succeed even without check_coherence
  fit <- suppressWarnings(fit_fewclust(
    d, y ~ day + arm, id = id, time = day, method = "mmm",
    control = fewclust_control(check_coherence = FALSE)
  ))
  expect_equal(fit$method, "mmm")
  expect_true(all(is.finite(tidy(fit)$p.value)))
  expect_gt(fit$convergence$trend$logLik, fit$convergence$constant$logLik - 50)
})

test_that("fc_mmm_negloglik guards against overflowing parameters without error", {
  d <- fc_test_data(15, seed = 2)
  des <- mmm_design(d, y ~ time + arm, trend = TRUE, transition = TRUE)
  val <- fewclust:::fc_mmm_negloglik(c(0, 0, 0, 0, 1, 500), des)
  expect_true(is.numeric(val) && length(val) == 1L)
})

test_that("Hessian inversion recovers finite beta SEs at the sigma_0 -> 0 boundary", {
  d <- simulate_clustered_binary(
    n_groups = 25, n_occasions = 4, effect = 0.4, sigma_u = 0, seed = 3
  )
  des <- mmm_design(d, y ~ time + arm, trend = FALSE, transition = FALSE)
  g <- glm(y ~ time + arm, family = binomial(), data = d)
  th <- c(coef(g) * des$scale, -12)
  V <- fewclust:::fc_mmm_vcov(th, des)
  se_beta <- sqrt(diag(V)[seq_len(des$p)])
  se_glm <- unname(summary(g)$coefficients[, "Std. Error"])
  expect_true(all(is.finite(se_beta)))
  expect_equal(unname(se_beta), se_glm, tolerance = 1e-2)
})
