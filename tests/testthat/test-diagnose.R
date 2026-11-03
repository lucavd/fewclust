test_that("diagnosis measures groups, occasions, prevalence and df", {
  dg <- diag_of(y ~ time * arm + x, fc_trial, group = id, time = time)
  expect_s3_class(dg, "cluster_diagnosis")
  expect_equal(dg$n_obs, 96L)
  expect_equal(dg$n_groups, 24L)
  expect_equal(dg$n_occasions, 4)
  expect_equal(dg$p_fixed, 5L)
  expect_equal(dg$df_t, 19L)
  expect_equal(dg$prevalence, mean(fc_trial$y))
  expect_equal(dg$separation$status, "none")
  expect_true(is.na(dg$icc) || (dg$icc > -1 && dg$icc < 1))
  expect_output(print(dg), "t reference on 19 df")
})

test_that("occasions are the longest series, with or without `time`", {
  d <- fc_test_data(10, seed = 1)
  d <- d[-c(1, 5, 9), ]
  dg <- diag_of(y ~ time + arm, d, group = id)
  expect_equal(dg$n_occasions, 4)
  short <- d[d$time < 1, ]
  short$time[short$id == "1"] <- short$time[short$id == "1"] + 10
  expect_equal(diag_of(y ~ arm, short, group = id, time = time)$n_occasions, 3)
  expect_equal(dg$event, "1")
  expect_null(dg$time)
  expect_output(print(dg), "no `time` given")
})

test_that("group and time accept bare names, strings and formulas", {
  g <- "id"
  d1 <- diag_of(y ~ time + arm, fc_trial, group = id, time = "time")
  d2 <- diag_of(y ~ time + arm, fc_trial, group = g, time = ~time)
  expect_equal(d1$n_groups, d2$n_groups)
  expect_equal(d1$n_occasions, d2$n_occasions)
  expect_error(diag_of(y ~ time + arm, fc_trial, group = nope),
               "not found")
  expect_error(diag_of(y ~ time + arm, fc_trial, group = id,
                                time = c("a", "b")), "exactly one column")
})

test_that("a random-effects term names the group and the slope structure", {
  dg <- diag_of(y ~ time + arm + (1 + time | id), fc_trial)
  expect_equal(dg$group, "id")
  expect_equal(dg$slope_vars, "time")
  expect_equal(dg$q_ranef, 2L)
  expect_true(dg$ranef_identifiable)
  two <- fc_test_data(20, seed = 1, n_occasions = 2)
  dg2 <- diag_of(y ~ time + arm + (1 + time | id), two)
  expect_false(dg2$ranef_identifiable)
  expect_error(diag_of(y ~ time + (1 | id) + (1 | arm), fc_trial),
               "single random-effects term")
  expect_error(diag_of(y ~ time + (1 | id), fc_trial, group = arm),
               "single grouping factor")
})

test_that("parsing validates the outcome, the design and missing data", {
  d <- fc_trial
  d$y[1:3] <- NA
  dg <- diag_of(y ~ time + arm, d, group = id)
  expect_equal(dg$n_dropped, 3L)
  expect_output(print(dg), "3 incomplete rows dropped")
  d$y <- d$y + 1
  expect_error(diag_of(y ~ time + arm, d, group = id), "binary")
  d <- fc_trial
  d$arm2 <- d$arm
  expect_error(diag_of(y ~ arm + arm2, d, group = id), "rank deficient")
  expect_error(diag_of(y ~ time + arm, fc_trial), "no grouping")
  d$time <- as.character(d$time)
  expect_error(diag_of(y ~ arm, d, group = id, time = time), "numeric")
})

test_that("logical and factor outcomes are accepted", {
  d <- fc_trial
  d$y <- d$y == 1
  expect_equal(diag_of(y ~ arm, d, group = id)$prevalence,
               mean(fc_trial$y))
  d$y <- factor(ifelse(fc_trial$y == 1, "yes", "no"), levels = c("no", "yes"))
  dg <- diag_of(y ~ arm, d, group = id)
  expect_equal(dg$prevalence, mean(fc_trial$y))
  expect_equal(dg$event, "yes")
  expect_output(print(dg), 'event "yes"')
})
