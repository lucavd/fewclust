test_that("control returns the study defaults", {
  ctrl <- fewclust_control()
  expect_s3_class(ctrl, "fewclust_control")
  expect_equal(ctrl$q, 20L)
  expect_true(ctrl$check)
  expect_true(ctrl$check_coherence)
  expect_equal(ctrl$check_threshold, 1e-3)
  expect_equal(ctrl$fg_bound, 0.75)
})

test_that("control validates its arguments", {
  expect_error(fewclust_control(q = 1), "at least 2")
  expect_error(fewclust_control(check = "yes"), "TRUE or FALSE")
  expect_error(fewclust_control(check_coherence = NA), "TRUE or FALSE")
  expect_error(fewclust_control(check_threshold = 0), "positive")
  expect_error(fewclust_control(fg_bound = 1), "strictly inside")
  expect_error(fewclust_control(corstr = "ar1"), "unused argument")
})
