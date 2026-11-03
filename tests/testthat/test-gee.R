test_that("fit_gee builds a broom-shaped coefficient table with t reference", {
  d <- fc_test_data(30, seed = 2)
  parsed <- fewclust:::fc_parse(y ~ time + arm + x, d,
                                group_expr = quote(id))
  res <- fewclust:::fit_gee(parsed, control = fewclust_control(),
                            type = "fg", reference = "t", alpha = 0.05,
                            df = 28L)
  expect_true(res$ok)
  tbl <- res$coefficients
  expect_named(tbl, c("term", "estimate", "std.error", "statistic", "df",
                      "p.value", "conf.low", "conf.high"))
  expect_equal(tbl$df, rep(28, 4))
  # t-based p and CI reproduce by hand
  expect_equal(tbl$p.value,
               2 * pt(-abs(tbl$statistic), df = 28))
  expect_equal(tbl$conf.high - tbl$estimate,
               qt(0.975, 28) * tbl$std.error)
})

test_that("z reference uses the normal distribution", {
  d <- fc_test_data(30, seed = 2)
  parsed <- fewclust:::fc_parse(y ~ time + arm, d, group_expr = quote(id))
  res <- fewclust:::fit_gee(parsed, control = fewclust_control(),
                            type = "unc", reference = "z", alpha = 0.05,
                            df = 28L)
  tbl <- res$coefficients
  expect_true(all(is.infinite(tbl$df)))
  expect_equal(tbl$p.value, 2 * pnorm(-abs(tbl$statistic)))
})

test_that("t and z references genuinely differ at small G", {
  d <- fc_test_data(12, seed = 6)
  parsed <- fewclust:::fc_parse(y ~ time + arm, d, group_expr = quote(id))
  rt <- fewclust:::fit_gee(parsed, fewclust_control(), "fg", "t", 0.05, 10L)
  rz <- fewclust:::fit_gee(parsed, fewclust_control(), "fg", "z", 0.05, 10L)
  expect_true(all(rt$coefficients$p.value >= rz$coefficients$p.value))
  expect_gt(min(rt$coefficients$p.value / rz$coefficients$p.value), 1)
})

test_that("a correction failure surfaces as a structured failure, not an error", {
  d <- fc_test_data(6, seed = 9)
  d$w <- as.integer(d$id == levels(d$id)[1L])
  d$y[d$id == levels(d$id)[1L]] <- c(0L, 1L, 0L, 1L)  # mixed: no separation
  parsed <- fewclust:::fc_parse(y ~ time + w, d, group_expr = quote(id))
  res <- fewclust:::fit_gee(parsed, fewclust_control(), "kc", "t", 0.05, 5L)
  expect_false(res$ok)
  expect_s3_class(res, "fc_failure")
  expect_match(res$reason, "leverage")
})

test_that("a 2-level factor outcome fits identically to its 0/1 recoding", {
  d <- fc_test_data(50, seed = 2)
  df <- d
  df$y <- factor(ifelse(d$y == 1, "yes", "no"), levels = c("no", "yes"))
  pn <- fewclust:::fc_parse(y ~ time + arm, d, group_expr = quote(id))
  pf <- fewclust:::fc_parse(y ~ time + arm, df, group_expr = quote(id))
  rn <- fewclust:::fit_gee(pn, fewclust_control(), "fg", "t", 0.05, 48L)
  rf <- fewclust:::fit_gee(pf, fewclust_control(), "fg", "t", 0.05, 48L)
  expect_true(rf$ok)
  expect_equal(rf$coefficients$estimate, rn$coefficients$estimate,
               tolerance = 1e-8)
})

test_that("fit warnings are captured into notes, not destroyed", {
  set.seed(41)
  d <- fc_test_data(20, seed = 41)
  d$wsep <- rnorm(nrow(d))
  d$y <- as.integer(d$wsep > 0)  # boundary fit: geeglm warns
  parsed <- fewclust:::fc_parse(y ~ wsep + arm, d, group_expr = quote(id))
  res <- fewclust:::fit_gee(parsed, fewclust_control(), "unc", "t", 0.05,
                            18L)
  if (isTRUE(res$ok)) {
    expect_match(paste(res$notes, collapse = " "), "GEE fit:")
  } else {
    expect_s3_class(res, "fc_failure")
  }
})

# The GEE strategy of the study (branch morris-v4, R/s7d/s7d_common.R::one_cor):
# geepack::geeglm(..., corstr = "ar1", scale.fix = TRUE) on data sorted by
# subject and time, with fewclust's Fay-Graubard variance.
test_that("the GEE strategy is the study's fit, whatever the order of the rows", {
  d <- fc_test_data(25, seed = 12, n_occasions = 5)
  d <- d[d$time < 1 | as.integer(as.character(d$id)) %% 3 != 0, ]   # some truncated series
  ref_data <- d[order(d$id, d$time), ]
  ref <- geepack::geeglm(y ~ time * arm, id = id, data = ref_data,
                         family = binomial, corstr = "ar1", scale.fix = TRUE)
  ref_se <- sqrt(diag(gee_vcov(ref, "fg")))
  shuffled <- d[sample(nrow(d)), ]
  fit <- fit_fewclust(shuffled, y ~ time * arm, id = id, time = time,
                      method = "gee")
  expect_equal(fit$method, "gee")
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(tidy(fit)$std.error, unname(ref_se), tolerance = 1e-10)
  expect_equal(unname(fit$engine$geese$alpha), unname(ref$geese$alpha), tolerance = 1e-10)
  expect_true(fit$engine$geese$model$scale.fix)
})
