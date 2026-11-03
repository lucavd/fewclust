# Reference implementation: the code of the study that produced the sign-flip
# results of the paper (branch morris-v4: R/s7c/s7c_tests.R::all_signs,
# R/s7d2/s7d2_ci.R::flip_signs and flip_p_b with diverge = TRUE, at b0 = 0),
# copied here with the data passed as plain vectors.  flip_test() must give
# the same p-value from the same RNG state.
study_flip_p <- function(id, y, x, X0, B = 1999L, enum_max = 12L) {
  id <- as.integer(factor(id))
  N <- length(unique(id))
  all_signs <- function(N) {
    m <- as.matrix(expand.grid(rep(list(c(1, -1)), N)))
    dimnames(m) <- NULL
    m
  }
  S <- if (N <= enum_max) list(G = all_signs(N), enum = TRUE) else
    list(G = matrix(sample(c(-1, 1), B * N, replace = TRUE), B, N), enum = FALSE)
  f <- tryCatch(suppressWarnings(glm.fit(X0, y, family = binomial(),
                                         control = glm.control(epsilon = 1e-12, maxit = 100))),
                error = function(e) NULL)
  ok <- !is.null(f) && isTRUE(f$converged) && all(is.finite(f$coefficients)) &&
    all(is.finite(f$fitted.values))
  if (!ok) return(list(p = NA_real_, why = "reduced fit failed"))
  if (max(abs(f$coefficients)) >= 1e3) return(list(p = NA_real_, why = "reduced fit diverged"))
  mu <- f$fitted.values; e <- y - mu; w <- mu * (1 - mu)
  I00 <- crossprod(X0, X0 * w); Ix0 <- colSums(X0 * (w * x))
  proj <- tryCatch(solve(I00, Ix0), error = function(e) NULL)
  if (is.null(proj)) return(list(p = NA_real_, why = "information singular"))
  u <- as.numeric(rowsum(x * e, id, reorder = TRUE)[, 1] -
                    rowsum(X0 * e, id, reorder = TRUE) %*% proj)
  if (!(sum(u^2) > 0)) return(list(p = NA_real_, why = "zero score variance"))
  if (max(abs(u)) < 1e-6) return(list(p = NA_real_, why = "score numerically zero"))
  To <- abs(sum(u)); Ts <- abs(as.numeric(S$G %*% u))
  p <- if (S$enum) mean(Ts >= To * (1 - 1e-10)) else (1 + sum(Ts >= To * (1 - 1e-10))) / (B + 1)
  list(p = p, why = NA_character_)
}

ref_on <- function(d, seed) {
  d <- d[order(d$id, d$time), ]
  set.seed(seed)
  study_flip_p(d$id, d$y, d$time * d$arm, cbind(1, d$time, d$arm))
}

test_that("flip_test reproduces the study's p-value with enumeration (N <= 12)", {
  for (s in 1:5) {
    d <- fc_test_data(10, seed = s, n_occasions = 4)
    ft <- flip_test(d, y ~ time * arm, id = id, term = "time:arm", seed = 100 + s)
    expect_s3_class(ft, "fewclust_flip")
    expect_true(ft$enumerated)
    expect_equal(ft$n_signs, 2^10)
    expect_identical(ft$p.value, ref_on(d, 100 + s)$p)
  }
})

test_that("flip_test reproduces the study's p-value with random signs (N > 12)", {
  for (s in 1:5) {
    d <- fc_test_data(25, seed = s, n_occasions = 3)
    ft <- flip_test(d, y ~ time * arm, id = id, term = "time:arm", seed = 200 + s)
    expect_false(ft$enumerated)
    expect_equal(ft$n_signs, 2000)
    expect_identical(ft$p.value, ref_on(d, 200 + s)$p)
  }
})

test_that("row order and column order of the data do not change the result", {
  d <- fc_test_data(25, seed = 7, n_occasions = 4)
  a <- flip_test(d, y ~ time * arm, id = id, time = time, term = "time:arm",
                 seed = 1)
  shuffled <- d[sample(nrow(d)), rev(names(d))]
  # with `time` the rows are put in the order of the study (group, then time)
  b <- flip_test(shuffled, y ~ time * arm, id = id, time = time,
                 term = "time:arm", seed = 1)
  expect_identical(a$p.value, b$p.value)
  b2 <- flip_test(shuffled, y ~ time * arm, id = id, term = "time:arm", seed = 1)
  expect_equal(b2$p.value, a$p.value)
  # the order of the terms in the formula does not matter either
  c2 <- flip_test(d, y ~ arm * time, id = id, term = "arm:time", seed = 1)
  expect_equal(c2$p.value, a$p.value)
})

test_that("the RNG state is used when no seed is given, and left alone when one is", {
  d <- fc_test_data(30, seed = 3)
  set.seed(42)
  a <- flip_test(d, y ~ time * arm, id = id, term = "time:arm")
  set.seed(42)
  b <- flip_test(d, y ~ time * arm, id = id, term = "time:arm")
  expect_identical(a$p.value, b$p.value)
  set.seed(9)
  before <- runif(1)
  set.seed(9)
  flip_test(d, y ~ time * arm, id = id, term = "time:arm", seed = 1)
  expect_identical(runif(1), before)
})

test_that("the default term is the last column of the mean model", {
  d <- fc_test_data(15, seed = 2)
  ft <- flip_test(d, y ~ time * arm, id = id, seed = 1)
  expect_equal(ft$term, "time:arm")
  expect_error(flip_test(d, y ~ time * arm, id = id, term = "nope"),
               "must be one of")
  expect_error(flip_test(d, y ~ 1, id = id), "no coefficient to test")
})

test_that("a covariate enters the model without the tested term", {
  d <- fc_test_data(20, seed = 4)
  ft <- flip_test(d, y ~ time * arm + x, id = id, term = "time:arm", seed = 3)
  d <- d[order(d$id, d$time), ]
  set.seed(3)
  ref <- study_flip_p(d$id, d$y, d$time * d$arm, cbind(1, d$time, d$arm, d$x))
  expect_identical(ft$p.value, ref$p)
})

test_that("an arm with only zero outcomes gives no answer (score numerically zero)", {
  d <- fc_test_data(20, seed = 6)
  d$y[d$arm == 1] <- 0L
  expect_warning(ft <- flip_test(d, y ~ time * arm, id = id, term = "time:arm", seed = 1), "not computable", class = "fewclust_warning")
  expect_false(ft$computable)
  expect_true(is.na(ft$p.value))
  expect_match(ft$reason, "diverged|numerically zero")
  expect_identical(ft$reason, ref_on(d, 1)$why)
  expect_output(print(ft), "not computable")
})

test_that("the test refuses an outcome that is not binary and needs two groups", {
  d <- fc_test_data(12, seed = 1)
  d$y <- d$y + 1L
  expect_error(flip_test(d, y ~ time * arm, id = id), "binary")
  one <- fc_test_data(12, seed = 1)
  one$g <- 1L
  expect_error(flip_test(one, y ~ time * arm, id = g), "fewer than 2 groups")
})

test_that("print says it is a p-value only", {
  ft <- flip_test(fc_trial, y ~ time * arm, id = id, seed = 1)
  expect_output(print(ft), "p-value only")
  expect_output(print(ft), "1,999 random sign assignments plus the identity")
  expect_true(ft$p.value > 0 && ft$p.value <= 1)
})

test_that("a tested term that is identically zero gives no answer, not an error", {
  d <- fc_test_data(12, seed = 3, n_occasions = 2)
  d <- d[!(d$arm == 1 & d$time > 0), ]   # no treated group after the first occasion
  expect_warning(ft <- flip_test(d, y ~ time * arm, id = id, time = time, term = "time:arm", seed = 1), "not computable", class = "fewclust_warning")
  expect_false(ft$computable)
  expect_identical(ft$reason, ref_on(d, 1)$why)
  expect_identical(ft$reason, "zero score variance")
  # the fitting function still refuses the rank-deficient design
  expect_error(suppressWarnings(fit_fewclust(d, y ~ time * arm, id = id, time = time)),
               "rank deficient", class = "fewclust_error")
})

# Decision D-82 of the study (branch morris-v4): a dataset with a constant
# outcome gives no usable result for any method.  The automatic strategies
# refuse it (rule D2); flip_test() says so instead of "score numerically zero".
test_that("a constant outcome is not computable, with its own reason", {
  for (val in 0:1) {
    d <- fc_test_data(12, seed = 4, n_occasions = 3)
    d$y <- val
    expect_warning(ft <- flip_test(d, y ~ time * arm, id = id, time = time, seed = 1), "not computable", class = "fewclust_warning")
    expect_true(is.na(ft$p.value))
    expect_false(ft$computable)
    expect_match(ft$reason, "^constant outcome")
    expect_output(print(ft), "not computable: constant outcome")
    expect_equal(ft$n_signs, 2^12)
    # the fitting function refuses the same dataset
    expect_error(fit_fewclust(d, y ~ time * arm, id = id, time = time),
                 "constant", class = "fewclust_error")
  }
  # the random stream is the same whether or not the outcome is constant
  d <- fc_test_data(30, seed = 4, n_occasions = 3)
  d0 <- d
  d0$y <- 0L
  set.seed(7)
  suppressWarnings(flip_test(d0, y ~ time * arm, id = id, time = time))
  after_const <- runif(1)
  set.seed(7)
  invisible(flip_test(d, y ~ time * arm, id = id, time = time))
  after_var <- runif(1)
  expect_identical(after_const, after_var)
  # one event is enough for the rule not to apply
  d1 <- d0
  d1$y[1L] <- 1L
  expect_false(isTRUE(startsWith(suppressWarnings(flip_test(d1, y ~ time * arm, id = id, time = time,
                                           seed = 1))$reason, "constant outcome")))
})

test_that("tidy() gives a one-row tibble", {
  ft <- flip_test(fc_trial, y ~ time * arm, id = id, time = time, seed = 1)
  td <- tidy(ft)
  expect_s3_class(td, "tbl_df")
  expect_named(td, c("term", "statistic", "p.value", "n_groups", "n_signs", "method"))
  expect_equal(td$p.value, ft$p.value)
})

test_that("the number of random sign assignments is set with `reps`", {
  d <- fc_test_data(20, seed = 2)
  ft <- flip_test(d, y ~ time * arm, id = id, reps = 499, seed = 1)
  expect_equal(ft$n_signs, 500)
})
