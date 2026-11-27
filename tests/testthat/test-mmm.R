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

test_that("the convergence check catches a false optimum and indefinite Hessian on real data", {
  d <- simulate_clustered_binary(n_groups = 10, n_occasions = 6, seed = 8)
  d$time <- round(d$time * 5)
  fit <- fit_fewclust(d, y ~ time * arm, id = id, time = time, method = "mmm")
  expect_equal(fit$method, "mmm")
  expect_true(fit$convergence$boundary)
  expect_equal(fit$convergence$latent_sd, "constant")
  expect_match(paste(fit$convergence$reason, collapse = " "), "false optimum")
  expect_match(
    paste(fit$convergence$reason, collapse = " "),
    "observed information matrix is not positive definite"
  )
  # With check_coherence = FALSE alone, check_hessian = TRUE still rejects it
  no_coh <- fit_fewclust(
    d, y ~ time * arm, id = id, time = time, method = "mmm",
    control = fewclust_control(check_coherence = FALSE)
  )
  expect_true(no_coh$convergence$boundary)
  expect_equal(no_coh$convergence$latent_sd, "constant")
  expect_match(
    paste(no_coh$convergence$reason, collapse = " "),
    "observed information matrix is not positive definite"
  )
  # Disabling both package extensions reproduces the published two-condition rule
  raw <- fit_fewclust(
    d, y ~ time * arm, id = id, time = time, method = "mmm",
    control = fewclust_control(check_coherence = FALSE, check_hessian = FALSE)
  )
  expect_false(raw$convergence$boundary)
  expect_equal(raw$convergence$latent_sd, "time trend")
  expect_lt(raw$engine$logLik, fit$engine$logLik)
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

test_that("fc_mmm_negloglik guards against overflowing parameters without error", {
  d <- fc_test_data(15, seed = 2)
  des <- mmm_design(d, y ~ time + arm, trend = TRUE, transition = TRUE)
  val <- fewclust:::fc_mmm_negloglik(c(0, 0, 0, 0, 1, 500), des)
  expect_true(is.numeric(val) && length(val) == 1L && is.infinite(val))
  expect_true(all(is.na(attr(val, "gradient"))))
})

test_that("fc_mmm_fit_once rejects a fit stopped on the a > 50 barrier", {
  d <- fc_test_data(15, seed = 2)
  des <- mmm_design(d, y ~ time + arm, trend = TRUE, transition = TRUE)
  local_mocked_bindings(
    nlm = function(f, p, ...) {
      list(
        minimum = 10, estimate = c(0, 0, 0, 0, 50, 1),
        gradient = numeric(6), code = 2L, iterations = 3L
      )
    },
    .package = "fewclust"
  )
  res <- fewclust:::fc_mmm_fit_once(des, numeric(6))
  expect_equal(res$code, 5L)
  expect_true(is.na(res$logLik))
  expect_true(all(is.na(res$vcov)))
})

test_that("long series scale the initial latent SD slope and converge cleanly", {
  d <- simulate_clustered_binary(n_groups = 20, n_occasions = 25, seed = 9)
  fit <- fit_fewclust(d, y ~ time + arm, id = id, time = time, method = "mmm")
  expect_equal(fit$method, "mmm")
  expect_true(fit$engine$code %in% c(1L, 2L))
  expect_true(all(is.finite(tidy(fit)$p.value)))
})

test_that("boundary block inversion at sigma_0 < 1e-3 recovers glm SEs and records a note", {
  d <- simulate_clustered_binary(
    n_groups = 25, n_occasions = 4, effect = 0.4, sigma_u = 0, seed = 3
  )
  des <- mmm_design(d, y ~ time + arm, trend = FALSE, transition = FALSE)
  g <- glm(y ~ time + arm, family = binomial(), data = d)
  th <- c(coef(g) * des$scale, -800)
  vc <- fewclust:::fc_mmm_vcov(th, des, check_hessian = TRUE)
  expect_true(vc$boundary_block)
  expect_false(vc$hessian_indefinite)
  se_beta <- sqrt(diag(vc$vcov)[seq_len(des$p)])
  se_glm <- unname(summary(g)$coefficients[, "Std. Error"])
  expect_equal(unname(se_beta), se_glm, tolerance = 1e-2)

  # Verify note is recorded on a fit where boundary_block activates
  parsed <- fewclust:::fc_parse(y ~ time + arm, d, quote(id))
  real_once <- fewclust:::fc_mmm_fit_once
  local_mocked_bindings(
    fc_mmm_fit_once = function(des_in, inits) {
      f <- real_once(des_in, inits)
      f$boundary_block <- TRUE
      f
    },
    .package = "fewclust"
  )
  fit <- suppressMessages(fit_fewclust(d, y ~ time + arm, id = id, method = "mmm"))
  expect_true(fit$convergence$boundary_block)
  expect_match(paste(fit$notes, collapse = " "), "conditional on log\\(sigma_0\\) held fixed")
})

test_that("continuous/jittered visit times scale cleanly and stay close to regular-grid fit", {
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
  expect_equal(f_jit$design$time_grid_regime, "continuous")
  expect_equal(f_jit$design$groups_with_gaps, 0L)
  expect_equal(f_jit$convergence$path, "trend fit kept")
  expect_true(f_jit$convergence$trend$code %in% c(1L, 2L))
  expect_lt(abs(f_jit$engine$logLik - f_reg$engine$logLik), 1.5)
  expect_lt(abs(unname(coef(f_jit)["arm"]) - unname(coef(f_reg)["arm"])), 0.1)
})

test_that("many unique jittered times converge cleanly without sub-model logLik collapse", {
  set.seed(1)
  d <- simulate_clustered_binary(n_groups = 120, n_occasions = 10, seed = 1)
  d$day <- d$time * 100 + runif(nrow(d), 0, 0.5)
  d <- d[order(d$id, d$day), ]
  fit <- suppressWarnings(fit_fewclust(
    d, y ~ day + arm, id = id, time = day, method = "mmm",
    control = fewclust_control(check_coherence = FALSE)
  ))
  expect_equal(fit$method, "mmm")
  expect_true(all(is.finite(tidy(fit)$p.value)))
  # Trend model nests constant model, so logLik must not collapse below constant
  expect_gt(fit$convergence$trend$logLik, fit$convergence$constant$logLik - 1e-2)
})

test_that("mmm_vcov = 'sandwich' computes exact Fay-Graubard sandwich SEs with and without time", {
  # With time
  fm <- suppressMessages(fit_fewclust(
    fc_trial, y ~ time * arm + x, id = id, time = time,
    method = "mmm", mmm_vcov = "model"
  ))
  rlang::reset_message_verbosity("fewclust_within_cluster")
  expect_no_message(
    fs <- fit_fewclust(
      fc_trial, y ~ time * arm + x, id = id, time = time,
      method = "mmm", mmm_vcov = "sandwich"
    )
  )
  expect_equal(fs$mmm_vcov, "sandwich")
  expect_equal(coef(fs), coef(fm), tolerance = 1e-10)
  expect_true(all(is.finite(tidy(fs)$std.error) & tidy(fs)$std.error > 0))
  expect_false(isTRUE(all.equal(tidy(fs)$std.error, tidy(fm)$std.error)))
  expect_output(print(fs), "Fay-Graubard sandwich")

  # Direct verification of B^{-1} sum_i (F_i U_i)(F_i U_i)' B^{-1} on constant-SD fit
  d_ord <- fc_trial[order(fc_trial$id), ]
  X <- model.matrix(y ~ arm + x, d_ord)
  des <- fewclust:::fc_mmm_design(
    d_ord$y, X, factor(d_ord$id), matrix(1, nrow(d_ord), 1),
    vcov_type = "sandwich", fg_bound = 0.75
  )
  fs0 <- suppressMessages(fit_fewclust(
    fc_trial, y ~ arm + x, id = id, method = "mmm", mmm_vcov = "sandwich"
  ))
  back <- c(1 / des$scale, 1)
  th <- unname(fs0$engine$theta / back)
  hc <- fewclust:::fc_mmm_hessian_clusters(th, des)
  Binv <- chol2inv(chol(hc$H))
  U_adj <- t(vapply(seq_len(nrow(hc$U)), function(i) {
    q <- diag(hc$H_list[[i]] %*% Binv)
    hc$U[i, ] / sqrt(1 - pmin(0.75, q))
  }, numeric(ncol(hc$U))))
  V_ref <- (Binv %*% crossprod(U_adj) %*% Binv) * tcrossprod(back)
  expect_equal(unname(vcov(fs0)), unname(V_ref[1:3, 1:3]), tolerance = 1e-8)

  # On a within-cluster random-slope design, sandwich SE exceeds model-based SE
  set.seed(77)
  G <- 40; m <- 10
  id_r <- factor(rep(seq_len(G), each = m))
  x_r <- rnorm(G * m)
  u0 <- rnorm(G, 0, 1)[as.integer(id_r)]
  u1 <- rnorm(G, 0, 2)[as.integer(id_r)]
  y_r <- rbinom(G * m, 1, plogis(-0.5 + u0 + u1 * x_r))
  dr <- data.frame(id = id_r, x = x_r, y = y_r)
  fm_r <- suppressMessages(fit_fewclust(dr, y ~ x, id = id, method = "mmm", mmm_vcov = "model"))
  fs_r <- suppressMessages(fit_fewclust(dr, y ~ x, id = id, method = "mmm", mmm_vcov = "sandwich"))
  expect_gt(tidy(fs_r)$std.error[2], 1.3 * tidy(fm_r)$std.error[2])
})
