test_that("generator returns the documented shape", {
  d <- simulate_clustered_binary(n_groups = 10, n_occasions = 3, seed = 1)
  expect_s3_class(d, "tbl_df")
  expect_named(d, c("id", "time", "arm", "x", "y"))
  expect_equal(nrow(d), 30L)
  expect_s3_class(d$id, "factor")
  expect_equal(nlevels(d$id), 10L)
  expect_true(all(d$y %in% 0:1))
  expect_setequal(unique(d$arm), c(0, 1))
})

test_that("same seed reproduces identical data, different seed does not", {
  a <- simulate_clustered_binary(seed = 42)
  b <- simulate_clustered_binary(seed = 42)
  expect_identical(a, b)
  d <- simulate_clustered_binary(seed = 43)
  expect_false(identical(a, d))
})

test_that("seed argument does not disturb the caller's RNG stream", {
  set.seed(7)
  before <- .Random.seed
  invisible(simulate_clustered_binary(seed = 1))
  expect_identical(.Random.seed, before)
})

test_that("random slopes vary across groups when slope_sd > 0", {
  slope_of <- function(dd) {
    vapply(split(dd, dd$id),
           function(g) coef(lm(y ~ time, data = g))[[2L]],
           numeric(1))
  }
  d1 <- simulate_clustered_binary(n_groups = 150, n_occasions = 40, slope_sd = 3,
                                  sigma_u = 0, effect = 0, seed = 5)
  d0 <- simulate_clustered_binary(n_groups = 150, n_occasions = 40, slope_sd = 0,
                                  sigma_u = 0, effect = 0, seed = 5)
  expect_gt(sd(slope_of(d1)), sd(slope_of(d0)))
})

test_that("prevalence is approximately honoured", {
  d <- simulate_clustered_binary(n_groups = 300, n_occasions = 4,
                                 prevalence = 0.3, effect = 0, seed = 9)
  expect_lt(abs(mean(d$y) - 0.3), 0.05)
})

test_that("invalid arguments are rejected", {
  expect_error(simulate_clustered_binary(n_groups = 1), "n_groups")
  expect_error(simulate_clustered_binary(prevalence = 0), "prevalence")
  expect_error(simulate_clustered_binary(n_occasions = 0), "n_occasions")
  expect_error(simulate_clustered_binary(sigma_u = -1), "sigma_u")
})
