test_that("print and summary describe the fit", {
  fit <- fit_fewclust(fc_trial, y ~ time + arm, id = id, time = time)
  expect_output(print(fit), "Marginalized multilevel model")
  expect_output(print(fit), "24 groups")
  s <- summary(fit)
  expect_s3_class(s, "summary.fewclust")
  expect_output(print(s), "Coefficients")
  expect_output(print(s), "t reference with 21 df")
  g <- fit_fewclust(fc_trial, y ~ time + arm, id = id, time = time, method = "gee")
  expect_output(print(g), "GEE")
  expect_output(print(g), "Fay-Graubard")
})

test_that("coef, vcov, confint and nobs follow base R conventions", {
  fit <- fit_fewclust(fc_trial, y ~ time + arm, id = id, time = time)
  expect_named(coef(fit), c("(Intercept)", "time", "arm"))
  expect_equal(dim(vcov(fit)), c(3L, 3L))
  expect_equal(unname(sqrt(diag(vcov(fit)))), tidy(fit)$std.error)
  ci <- confint(fit)
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(unname(ci[, 2] - coef(fit)), qt(0.975, fit$df) * tidy(fit)$std.error)
  expect_equal(rownames(confint(fit, "arm", level = 0.9)), "arm")
  expect_equal(colnames(confint(fit, level = 0.9)), c("5 %", "95 %"))
  expect_equal(nobs(fit), 96L)
  expect_equal(formula(fit), y ~ time + arm)
})

test_that("predict gives the marginal log-odds and probabilities", {
  fit <- fit_fewclust(fc_trial, y ~ time + arm, id = id, time = time)
  eta <- predict(fit)
  expect_length(eta, 96L)
  expect_equal(unname(eta), unname(drop(model.matrix(~ time + arm, fc_trial) %*% coef(fit))))
  expect_equal(predict(fit, type = "response"), plogis(eta))
  nd <- data.frame(time = c(0, 1), arm = c(1, 1))
  expect_equal(unname(predict(fit, newdata = nd)),
               drop(cbind(1, nd$time, nd$arm) %*% coef(fit)))
})

test_that("factor covariates are handled in predict", {
  d <- fc_test_data(30, seed = 3)
  d$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  fit <- fit_fewclust(d, y ~ time + grp, id = id, time = time)
  nd <- data.frame(time = 0.5, grp = factor("c", levels = c("a", "b", "c")))
  expect_equal(unname(predict(fit, newdata = nd)),
               unname(coef(fit)[1] + 0.5 * coef(fit)[2] + coef(fit)["grpc"]))
})
