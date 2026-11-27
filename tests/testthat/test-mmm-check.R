# The convergence check of the marginalized model on the six study datasets
# of fixtures/vedovelli (scenario M3, GenBinaryY, 25 clusters, 8 occasions,
# null arm; regenerated from their seeds on branch morris-v4 by
# R/s6b/s6b_vedovelli_fixtures.R).  expected.csv holds the study fit with the
# time trend of the latent standard deviation (MMLB::mm, fit A), the branch of
# the rule, the estimate and standard error used, and the p-value with a t
# reference on G - p = 25 - 4 = 21 degrees of freedom.  The rule is the one of
# the paper: threshold 1e-3, nlm codes other than 1 and 2, fallback to the
# constant latent standard deviation.

vedovelli <- function(file, control) {
  d <- read.csv(test_path("fixtures", "vedovelli", file))
  parsed <- fewclust:::fc_parse(binY ~ time * trt, d, quote(id), quote(time))
  fewclust:::fit_mmm(parsed, control, "t", 0.05, 21L)
}

test_that("the check reproduces the branch, estimate and p-value of the study", {
  ex <- read.csv(test_path("fixtures", "vedovelli", "expected.csv"))
  published <- fewclust_control(check_coherence = FALSE)
  for (k in seq_len(nrow(ex))) {
    e <- ex[k, ]
    res <- vedovelli(e$file, published)
    info <- e$file
    ck <- res$model$check
    tab <- res$coefficients
    j <- tab$term == "time:trt"
    expect_equal(ck$boundary, e$at_boundary, info = info)
    expect_equal(ck$path, if (e$used == "B") "constant fit used" else "trend fit kept",
                 info = info)
    expect_equal(ck$latent_sd, if (e$used == "B") "constant" else "time trend",
                 info = info)
    # fit A as stored by the study.  Where nlm did not converge (codes 3-5)
    # the point where it stops depends on the implementation (case_04, code
    # 3: estimate 0.0054328 here against 0.0054321 in MMLB, standard error
    # 0.086 against 0.099, both collapsed against 0.225 without the trend);
    # only the code is compared there.
    expect_equal(ck$trend$code, e$A_code, info = info)
    if (e$A_code %in% c(1L, 2L)) {
      expect_equal(unname(ck$trend$estimate["time:trt"]), e$A_est, tolerance = 1e-6,
                   info = info)
      expect_equal(unname(ck$trend$std.error["time:trt"]), e$A_se, tolerance = 1e-5,
                   info = info)
    }
    if (e$A_latent_sd_slope >= 1e-3) {
      expect_equal(ck$trend$slope_min, e$A_latent_sd_slope, tolerance = 1e-4,
                   info = info)
    } else {
      expect_lt(ck$trend$slope_min, 1e-3)
    }
    # the fit that is reported
    expect_equal(tab$estimate[j], e$est, tolerance = 1e-6, info = info)
    expect_equal(tab$std.error[j], e$se, tolerance = 1e-5, info = info)
    expect_equal(tab$df[j], e$df_t, info = info)
    expect_equal(tab$p.value[j], e$p_t, tolerance = 1e-5, info = info)
    # the reason is recorded when the check intervenes
    if (e$at_boundary) {
      expect_match(paste(ck$reason, collapse = " "),
                   if (grepl("nonconvergence", e$branch)) "did not converge"
                   else "slope of the latent standard deviation", info = info)
      expect_false(is.null(ck$constant))
    } else {
      expect_length(ck$reason, 0L)
    }
  }
})

test_that("the default check gives the same branch on the six study datasets", {
  ex <- read.csv(test_path("fixtures", "vedovelli", "expected.csv"))
  for (k in seq_len(nrow(ex))) {
    res <- vedovelli(ex$file[k], fewclust_control())
    expect_equal(res$model$check$boundary, ex$at_boundary[k], info = ex$file[k])
  }
})

# The branches of the rule of the study (D-90): the fit without the trend
# replaces a boundary fit only if nlm reports its convergence (code 1 or 2)
# and its estimates and standard errors are finite; otherwise no result.
# The refit is forced with a stub of the engine on the first study dataset
# (case_01: trend fit at the boundary, slope 1.2e-8, refit with code 1).

# the engine as it is, captured before any stub; a stub of
# fc_mmm_fit_once() alters the refit without the trend (d$m_s == 1) and
# leaves the trend fit as it is
real_fit_once <- fewclust:::fc_mmm_fit_once
stub_refit <- function(alter) {
  function(d, inits) {
    f <- real_fit_once(d, inits)
    if (d$m_s == 1L) alter(f) else f
  }
}

vedovelli_data <- function(k = 1L) {
  ex <- read.csv(test_path("fixtures", "vedovelli", "expected.csv"))
  read.csv(test_path("fixtures", "vedovelli", ex$file[k]))
}

study_rule <- fewclust_control(check_coherence = FALSE)

test_that("a converged refit without the trend is used", {
  ex <- read.csv(test_path("fixtures", "vedovelli", "expected.csv"))
  for (k in which(ex$used == "B")) {
    res <- vedovelli(ex$file[k], study_rule)
    ck <- res$model$check
    expect_true(res$ok, info = ex$file[k])
    expect_equal(ck$path, "constant fit used", info = ex$file[k])
    expect_equal(ck$refit_status, "converged", info = ex$file[k])
    expect_true(ck$constant$code %in% c(1L, 2L), info = ex$file[k])
    expect_true(res$model$code %in% c(1L, 2L), info = ex$file[k])
  }
  fit <- fit_fewclust(vedovelli_data(), binY ~ time * trt,
    id = id, time = time, method = "mmm", control = study_rule
  )
  expect_equal(fit$method, "mmm")
  expect_equal(fit$convergence$path, "constant fit used")
  expect_output(print(summary(fit)), "intervened \\(trend fit at the boundary\\)")
})

test_that("a refit with nlm code 3, 4 or 5 gives no result", {
  for (code in 3:5) {
    local_mocked_bindings(
      fc_mmm_fit_once = stub_refit(function(f) {
        f$code <- code
        f
      }),
      .package = "fewclust"
    )
    res <- vedovelli(read.csv(test_path(
      "fixtures", "vedovelli", "expected.csv"
    ))$file[1L], study_rule)
    expect_s3_class(res, "fc_failure")
    expect_equal(res$check$path, "no result")
    expect_equal(res$check$latent_sd, "no result")
    expect_equal(res$check$refit_status, "not converged")
    expect_true(res$check$boundary)
    expect_equal(res$check$constant$code, code)
    expect_match(res$reason, sprintf("did not converge \\(nlm code %d\\)", code))
    expect_match(res$reason, "slope of the latent standard deviation")
    expect_match(res$reason, "no result")
    # method = "mmm": an error of the class of every failure of the model
    err <- expect_error(
      fit_fewclust(vedovelli_data(), binY ~ time * trt,
        id = id, time = time, method = "mmm", control = study_rule
      ),
      sprintf("did not converge \\(nlm code %d\\)", code),
      class = "fewclust_error"
    )
    expect_equal(err$convergence$path, "no result")
  }
})

test_that("a trend fit that did not converge is not reported either", {
  # trend fit with code 4, refit with code 3: neither is reported
  local_mocked_bindings(
    fc_mmm_fit_once = function(d, inits) {
      f <- real_fit_once(d, inits)
      f$code <- if (d$m_s == 1L) 3L else 4L
      f
    },
    .package = "fewclust"
  )
  ex <- read.csv(test_path("fixtures", "vedovelli", "expected.csv"))
  res <- vedovelli(ex$file[5L], study_rule) # case_05: trend fit kept
  expect_s3_class(res, "fc_failure")
  expect_match(res$reason, "the optimiser did not converge \\(nlm code 4\\)")
  expect_match(res$reason, "did not converge \\(nlm code 3\\)")
})

test_that("a refit that is not finite gives no result", {
  ex <- read.csv(test_path("fixtures", "vedovelli", "expected.csv"))
  # a non-finite estimate, with a converged code
  local_mocked_bindings(
    fc_mmm_fit_once = stub_refit(function(f) {
      f$theta[2L] <- NA_real_
      f
    }),
    .package = "fewclust"
  )
  res <- vedovelli(ex$file[1L], study_rule)
  expect_s3_class(res, "fc_failure")
  expect_equal(res$check$path, "no result")
  expect_equal(res$check$refit_status, "not finite")
  expect_match(res$reason, "no finite estimate or standard error \\(nlm code 1\\)")
  # a non-finite standard error
  local_mocked_bindings(
    fc_mmm_fit_once = stub_refit(function(f) {
      f$vcov[] <- NA_real_
      f
    }),
    .package = "fewclust"
  )
  res <- vedovelli(ex$file[1L], study_rule)
  expect_s3_class(res, "fc_failure")
  expect_equal(res$check$refit_status, "not finite")
  # an error of the refit
  local_mocked_bindings(
    fc_mmm_fit_once = function(d, inits) {
      if (d$m_s == 1L) stop("stubbed failure of the constant fit")
      real_fit_once(d, inits)
    },
    .package = "fewclust"
  )
  res <- vedovelli(ex$file[1L], study_rule)
  expect_s3_class(res, "fc_failure")
  expect_equal(res$check$path, "no result")
  expect_equal(res$check$refit_status, "not finite")
  expect_null(res$check$constant)
  expect_match(res$reason, "stubbed failure of the constant fit")
})

test_that("with method = \"auto\" no result falls back to GEE", {
  local_mocked_bindings(
    fc_mmm_fit_once = stub_refit(function(f) {
      f$code <- 4L
      f
    }),
    .package = "fewclust"
  )
  d <- vedovelli_data()
  expect_message(
    fit <- fit_fewclust(d, binY ~ time * trt,
      id = id, time = time, control = study_rule
    ),
    "did not converge \\(nlm code 4\\)",
    class = "fewclust_message"
  )
  gee <- fit_fewclust(d, binY ~ time * trt,
    id = id, time = time, method = "gee"
  )
  expect_equal(fit$method, "gee")
  expect_match(fit$fallback, "nlm code 4")
  expect_identical(tidy(fit), tidy(gee))
  expect_identical(vcov(fit), vcov(gee))
  expect_equal(fit$convergence$path, "no result")
  expect_equal(fit$convergence$refit_status, "not converged")
  expect_null(gee$convergence)
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("marginalized multilevel model: no result", out)))
  expect_true(any(grepl("did not converge \\(nlm code 4\\)", out)))
})

test_that("a trend fit error still tries the constant-SD refit", {
  d <- vedovelli_data()
  local_mocked_bindings(
    fc_mmm_fit_once = function(des, inits) {
      if (des$m_s > 1L) stop("trend numerical error")
      real_fit_once(des, inits)
    },
    .package = "fewclust"
  )
  fit <- fit_fewclust(d, binY ~ time * trt, id = id, time = time, method = "mmm")
  expect_equal(fit$method, "mmm")
  expect_equal(fit$convergence$path, "constant fit used")
  expect_true(fit$convergence$boundary)
  expect_match(paste(fit$convergence$reason, collapse = " "), "trend numerical error")
})

test_that("a trend fit with code 1 and non-finite SE falls back to converged constant refit", {
  d <- vedovelli_data(5L)
  local_mocked_bindings(
    fc_mmm_fit_once = function(des, inits) {
      f <- real_fit_once(des, inits)
      if (des$m_s > 1L) {
        f$code <- 1L
        f$slope_min <- 0.5
        f$vcov[] <- NA_real_
      }
      f
    },
    .package = "fewclust"
  )
  fit <- fit_fewclust(d, binY ~ time * trt, id = id, time = time, method = "mmm")
  expect_equal(fit$method, "mmm")
  expect_equal(fit$convergence$path, "constant fit used")
  expect_true(fit$convergence$boundary)
  expect_true(all(is.finite(tidy(fit)$p.value)))
})

test_that("without time, a non-converged MMM is rejected when check = TRUE", {
  d <- fc_test_data(20, seed = 4)
  parsed <- fewclust:::fc_parse(y ~ time + arm, d, quote(id))
  local_mocked_bindings(
    fc_mmm_fit_once = function(des, inits) {
      f <- real_fit_once(des, inits)
      f$code <- 5L
      f
    },
    .package = "fewclust"
  )
  res <- fewclust:::fit_mmm(parsed, fewclust_control(check = TRUE), "t", 0.05, 17L)
  expect_s3_class(res, "fc_failure")
  expect_match(res$reason, "nlm code 5")
  # With check = FALSE the unconverged fit is still returned
  raw <- fewclust:::fit_mmm(parsed, fewclust_control(check = FALSE), "t", 0.05, 17L)
  expect_true(raw$ok)
  expect_equal(raw$model$code, 5L)
})
