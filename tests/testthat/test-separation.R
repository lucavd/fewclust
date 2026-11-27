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

test_that("quasi-complete separation in a 2x2 interaction cell is detected", {
  # Zero events in control arm at follow-up (time == 1, arm == 0):
  # separating vector is time - time:arm, so no single column separates alone
  d <- simulate_clustered_binary(
    n_groups = 20, n_occasions = 2, prevalence = 0.3, effect = 0, seed = 10
  )
  d$y[d$time == 1 & d$arm == 0] <- 0L
  dg <- diag_of(y ~ time * arm, d, group = id, time = time)
  expect_equal(dg$separation$status, "quasi")
  expect_error(
    fit_fewclust(d, y ~ time * arm, id = id, time = time),
    "separat",
    class = "fewclust_error"
  )
})

test_that("multi-column quasi-complete separation in an additive model is detected", {
  set.seed(11)
  n <- 80
  d <- data.frame(
    id = factor(rep(1:20, each = 4)),
    a = rep(c(0L, 1L), each = 2, times = 20),
    b = rep(c(0L, 1L), times = 40)
  )
  # Separating vector: -1 + a + b (<0 at (0,0), >0 at (1,1), 0 at (0,1)&(1,0))
  d$y <- ifelse(d$a == 0 & d$b == 0, 0L,
    ifelse(d$a == 1 & d$b == 1, 1L, rbinom(n, 1, 0.5))
  )
  dg <- diag_of(y ~ a + b, d, group = id)
  expect_equal(dg$separation$status, "quasi")
})

test_that("a single pure cell in an additive 2x2 model is NOT separation", {
  set.seed(12)
  n <- 120
  d <- data.frame(
    id = factor(rep(1:30, each = 4)),
    a = rep(c(0L, 1L), each = 2, times = 30),
    b = rep(c(0L, 1L), times = 60)
  )
  # Only (1,1) is pure; (0,0), (0,1), (1,0) are mixed -> finite MLE
  d$y <- ifelse(d$a == 1 & d$b == 1, 1L, rbinom(n, 1, 0.45))
  dg <- diag_of(y ~ a + b, d, group = id)
  expect_equal(dg$separation$status, "none")
})

test_that("separation detection is invariant to covariate scaling", {
  d <- simulate_clustered_binary(
    n_groups = 20, n_occasions = 2, prevalence = 0.3, effect = 0, seed = 10
  )
  d$y[d$time == 1 & d$arm == 0] <- 0L
  d$time_small <- d$time * 1e-5
  d$time_big <- d$time * 1e5
  expect_equal(diag_of(y ~ time_small * arm, d, group = id)$separation$status, "quasi")
  expect_equal(diag_of(y ~ time_big * arm, d, group = id)$separation$status, "quasi")

  d2 <- fc_trial
  d2$x_small <- d2$x * 1e-5
  d2$x_big <- d2$x * 1e5
  expect_equal(diag_of(y ~ time * arm + x_small, d2, group = id)$separation$status, "none")
  expect_equal(diag_of(y ~ time * arm + x_big, d2, group = id)$separation$status, "none")
})
