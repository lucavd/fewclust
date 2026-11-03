test_that("complete separation by a continuous covariate is detected pre-fit", {
  set.seed(21)
  d <- simulate_clustered_binary(n_groups = 20, n_occasions = 4, seed = 21)
  d$w <- rnorm(nrow(d))
  d$y <- as.integer(d$w > 0)
  dg <- diag_of(y ~ w + arm, d, group = id)
  expect_equal(dg$separation$status, "complete")
})

test_that("quasi-separation by a pure dummy cell is detected and named", {
  d <- simulate_clustered_binary(n_groups = 20, n_occasions = 4, seed = 22)
  d$y[d$arm == 1] <- 1L  # treated groups all respond
  dg <- diag_of(y ~ time + arm, d, group = id)
  expect_true(dg$separation$status %in% c("quasi", "complete"))
  expect_match(paste(dg$separation$evidence, collapse = " "), "arm")
})

test_that("well-behaved data are not flagged", {
  dg <- diag_of(y ~ time * arm + x, fc_trial, group = id)
  expect_equal(dg$separation$status, "none")
})

test_that("monotone separation by a discrete covariate is detected exactly", {
  d <- simulate_clustered_binary(n_groups = 30, n_occasions = 4, seed = 6)
  d$dose <- rep(0:2, length.out = nrow(d))
  # complete: dose >= 2 -> event, dose < 2 -> no event
  d$y <- as.integer(d$dose >= 2)
  dg <- diag_of(y ~ dose + arm, d, group = id)
  expect_equal(dg$separation$status, "complete")
  expect_match(paste(dg$separation$evidence, collapse = " "), "dose")
})

test_that("a pure non-extreme cell is NOT separation (the MLE is finite)", {
  set.seed(5)
  n <- 120
  d <- data.frame(
    id = factor(rep(1:30, each = 4)),
    dose = rep(0:2, length.out = n),
    x = rnorm(n)
  )
  # all events at the top dose, mixed outcomes below: glm converges to a
  # finite coefficient, so this must not be refused
  d$y <- ifelse(d$dose == 2, 1L, rbinom(n, 1, 0.45))
  dg <- diag_of(y ~ dose + x, d, group = id)
  expect_equal(dg$separation$status, "none")
})

test_that("the max_tab_levels threshold governs the boundary-tie check", {
  d <- simulate_clustered_binary(n_groups = 30, n_occasions = 4, seed = 6)
  d$dose <- rep(0:2, length.out = nrow(d))
  # quasi-complete: y = 1 iff dose == 2, except one dose-2 row flipped to
  # 0, so the outcome classes overlap only at the boundary value 2
  d$y <- as.integer(d$dose == 2)
  d$y[which(d$dose == 2)[1]] <- 0L
  X <- cbind("(Intercept)" = 1, dose = d$dose)
  default_dg <- fewclust:::detect_separation(X, d$y, max_tab_levels = 12L)
  narrow_dg <- fewclust:::detect_separation(X, d$y, max_tab_levels = 2L)
  tab_default <- grepl("boundary value", paste(default_dg$evidence,
                                               collapse = " "))
  tab_narrow <- grepl("boundary value", paste(narrow_dg$evidence,
                                              collapse = " "))
  expect_true(tab_default)
  expect_false(tab_narrow)
})

test_that("intercept-free designs are not falsely refused for separation", {
  set.seed(8)
  d <- fc_test_data(20, seed = 8)
  # ranges disjoint but strictly positive: with y ~ 0 + ... the MLE is
  # finite, so this must not be flagged
  d$age <- ifelse(d$y == 1, runif(nrow(d), 1, 2), runif(nrow(d), 3, 4))
  dg <- diag_of(y ~ 0 + age + time, d, group = id)
  expect_equal(dg$separation$status, "none")
  fit <- fit_fewclust(d, y ~ 0 + age + time, id = id, time = time)
  expect_s3_class(fit, "fewclust")
  # genuine through-origin separation IS still caught
  d$sgn <- ifelse(d$y == 1, abs(d$x) + 0.1, -abs(d$x) - 0.1)
  dg2 <- diag_of(y ~ 0 + sgn + time, d, group = id)
  expect_equal(dg2$separation$status, "complete")
})

test_that("a constant outcome yields no spurious separation evidence", {
  d <- simulate_clustered_binary(n_groups = 12, seed = 1)
  d$y <- 1L
  parsed_X <- cbind(1, time = d$time, arm = d$arm)
  res <- fewclust:::detect_separation(parsed_X, d$y)
  expect_equal(res$status, "none")
  expect_length(res$evidence, 0L)
})

test_that("near-boundary but non-separated data are not flagged", {
  # rare outcome, but no covariate pattern predicts perfectly
  d <- simulate_clustered_binary(n_groups = 60, n_occasions = 6,
                                 prevalence = 0.08, effect = 0.3, seed = 23)
  dg <- diag_of(y ~ time + arm + x, d, group = id)
  expect_equal(dg$separation$status, "none")
})
