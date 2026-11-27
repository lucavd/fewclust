test_that("the default fit is the marginalized model with a t reference", {
  fit <- fit_fewclust(fc_trial, y ~ time * arm + x, id = id, time = time)
  expect_s3_class(fit, "fewclust")
  expect_equal(fit$method, "mmm")
  expect_null(fit$fallback)
  expect_equal(fit$df, 24 - 5)
  tab <- tidy(fit)
  expect_equal(tab$df, rep(19, 5))
  expect_equal(tab$p.value, 2 * pt(-abs(tab$statistic), 19))
  expect_s3_class(fit$engine, "fewclust_mmm")
  expect_equal(fit$n_groups, 24L)
  expect_equal(fit$n_obs, 96L)
  expect_equal(fit$n_occasions, 4)
})

test_that("method = 'gee' is GEE with AR(1) and Fay-Graubard by default", {
  fit <- fit_fewclust(fc_trial, y ~ time * arm + x, id = id, time = time, method = "gee")
  expect_equal(fit$method, "gee")
  expect_equal(fit$corstr, "ar1")
  expect_equal(fit$correction, "fg")
  expect_s3_class(fit$engine, "geeglm")
  expect_equal(tidy(fit)$std.error, unname(sqrt(diag(gee_vcov(fit$engine, "fg")))))
})

test_that("corstr and correction choose the GEE variant", {
  d <- fc_test_data(30, seed = 2)
  ex <- fit_fewclust(d, y ~ time + arm, id = id, method = "gee", corstr = "exchangeable")
  expect_equal(ex$corstr, "exchangeable")
  kc <- suppressWarnings(fit_fewclust(d, y ~ time + arm, id = id, method = "gee", correction = "kc"))
  expect_equal(tidy(kc)$std.error, unname(sqrt(diag(gee_vcov(kc$engine, "kc")))))
  un <- suppressWarnings(fit_fewclust(d, y ~ time + arm, id = id, method = "gee", correction = "none"))
  expect_equal(tidy(un)$std.error, unname(sqrt(diag(gee_vcov(un$engine, "unc")))))
  expect_error(fit_fewclust(d, y ~ time + arm, id = id, method = "glmm"), "must be one of")
})

test_that("the default falls back to GEE when the marginalized model fails", {
  d <- fc_test_data(30, seed = 2)
  local_mocked_bindings(fit_mmm = function(...) fewclust:::fc_failure("stubbed failure"),
                        .package = "fewclust")
  expect_message(fit <- fit_fewclust(d, y ~ time + arm, id = id, time = time), "stubbed failure",
                 class = "fewclust_message")
  expect_equal(fit$method, "gee")
  expect_equal(fit$corstr, "ar1")
  expect_match(fit$fallback, "stubbed failure")
  # an explicit request does not fall back
  expect_error(fit_fewclust(d, y ~ time + arm, id = id, time = time, method = "mmm"), "stubbed failure",
               class = "fewclust_error")
})

test_that("when both models fail the error points to flip_test()", {
  d <- fc_test_data(30, seed = 2)
  local_mocked_bindings(fit_mmm = function(...) fewclust:::fc_failure("stub MMM"),
                        fc_gee_model = function(...) fewclust:::fc_failure("stub GEE"),
                        .package = "fewclust")
  err <- expect_error(fit_fewclust(d, y ~ time + arm, id = id, time = time), class = "fewclust_error")
  expect_match(conditionMessage(err), "stub MMM")
  expect_match(conditionMessage(err), "stub GEE")
  expect_match(conditionMessage(err), "flip_test")
})

test_that("without time, auto defaults to exchangeable GEE-FG and ar1 warns", {
  d <- fc_test_data(30, seed = 2)
  expect_message(
    fit <- fit_fewclust(d, y ~ arm + x, id = id),
    "exchangeable",
    class = "fewclust_message"
  )
  expect_equal(fit$method, "gee")
  expect_equal(fit$corstr, "exchangeable")
  expect_equal(fit$correction, "fg")

  g <- suppressMessages(fit_fewclust(d, y ~ arm + x, id = id, method = "gee"))
  expect_equal(g$corstr, "exchangeable")

  expect_warning(
    g_ar1 <- suppressMessages(
      fit_fewclust(d, y ~ arm + x, id = id, method = "gee", corstr = "ar1")
    ),
    "order of",
    class = "fewclust_warning"
  )
  expect_equal(g_ar1$corstr, "ar1")

  # explicit method = "mmm" without time still fits the constant-SD MMM
  m <- suppressMessages(fit_fewclust(d, y ~ arm + x, id = id, method = "mmm"))
  expect_equal(m$method, "mmm")
  expect_equal(m$convergence$latent_sd, "constant")
})

test_that("inference is invariant to covariate unit conversions", {
  d <- fc_trial
  d$time_yr <- d$time / 50
  d$x_small <- d$x / 100
  d$x_big <- d$x * 1000
  f0 <- fit_fewclust(d, y ~ time * arm + x, id = id, time = time)
  f_yr <- fit_fewclust(d, y ~ time_yr * arm + x_small, id = id, time = time_yr)
  f_big <- fit_fewclust(d, y ~ time * arm + x_big, id = id, time = time)
  t0 <- tidy(f0)
  ty <- tidy(f_yr)
  tb <- tidy(f_big)
  expect_true(all(is.finite(ty$p.value)))
  expect_true(all(is.finite(tb$p.value)))
  expect_equal(ty$statistic, t0$statistic, tolerance = 1e-3)
  expect_equal(ty$p.value, t0$p.value, tolerance = 1e-3)
  expect_equal(tb$statistic, t0$statistic, tolerance = 1e-3)
  expect_equal(tb$p.value, t0$p.value, tolerance = 1e-3)
})

test_that("a coefficient failing the standardized usability criterion gets NA inference and a warning", {
  d <- fc_test_data(40, seed = 3, n_occasions = 6)
  d$x_big <- d$x * 1000
  real_fit <- fewclust:::fc_mmm_fit_once
  local_mocked_bindings(
    fc_mmm_fit_once = function(des, inits) {
      f <- real_fit(des, inits)
      # Inflate the SE of x_big to 12 on the standardized scale (0.012 raw)
      k <- match("x_big", names(f$theta))
      s <- sd(d$x_big)
      f$vcov[k, k] <- (12 / s)^2
      f
    },
    .package = "fewclust"
  )
  expect_warning(
    fit <- fit_fewclust(d, y ~ time * arm + x_big, id = id, time = time, method = "mmm"),
    "x_big",
    class = "fewclust_warning"
  )
  tab <- tidy(fit)
  bad <- tab$term == "x_big"
  expect_true(is.na(tab$p.value[bad]))
  expect_true(all(is.finite(tab$p.value[!bad])))
})

test_that("incomplete rows are dropped and counted", {
  d <- fc_trial
  d$x[1:3] <- NA
  fit <- fit_fewclust(d, y ~ time + arm + x, id = id, time = time)
  expect_equal(fit$n_obs, 93L)
  expect_match(paste(fit$notes, collapse = " "), "3 rows with missing values")
})

test_that("conf_level sets the intervals", {
  fit <- fit_fewclust(fc_trial, y ~ time + arm, id = id, time = time, conf_level = 0.9)
  tab <- tidy(fit, conf.int = TRUE, conf.level = 0.9)
  expect_equal(tab$conf.high - tab$estimate, qt(0.95, fit$df) * tab$std.error)
  expect_error(fit_fewclust(fc_trial, y ~ arm, id = id, conf_level = 1), "conf_level")
})

test_that("argument validation gives clear errors", {
  expect_error(fit_fewclust(fc_trial, y ~ arm), "id", class = "fewclust_error")
  expect_error(fit_fewclust(as.matrix(fc_trial), y ~ arm, id = id), "data frame",
               class = "fewclust_error")
  expect_error(fit_fewclust(fc_trial, ~ arm, id = id), "two-sided", class = "fewclust_error")
  expect_error(fit_fewclust(fc_trial, y ~ arm, id = nope))
  expect_error(fit_fewclust(fc_trial, y ~ arm, id = id, control = list()), "fewclust_control")
  expect_error(fit_fewclust(fc_trial, y ~ arm + (1 | id), id = id), "random-effects",
               class = "fewclust_error")
})
