test_that("tibbles, pipes and bare or quoted column names work", {
  tb <- tibble::as_tibble(fc_trial)
  a <- tb |> fit_fewclust(y ~ time * arm, id = id, time = time)
  b <- fit_fewclust(fc_trial, y ~ time * arm, id = "id", time = "time")
  expect_equal(coef(a), coef(b))
})

test_that("column names can be passed through a user function with {{ }}", {
  my_fit <- function(data, cluster, occasion) {
    fit_fewclust(data, y ~ time * arm, id = {{ cluster }}, time = {{ occasion }})
  }
  a <- my_fit(fc_trial, id, time)
  b <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
  expect_equal(coef(a), coef(b))
  expect_error(fit_fewclust(fc_trial, y ~ time, id = c(id, arm)), "one column")
})

test_that("tidy, glance and augment return tibbles with broom's names", {
  fit <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_named(td, c("term", "estimate", "std.error", "statistic", "df", "p.value"))
  expect_named(tidy(fit, conf.int = TRUE),
               c("term", "estimate", "std.error", "statistic", "df", "p.value",
                 "conf.low", "conf.high"))
  gl <- glance(fit)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("nobs", "n_groups", "method", "df") %in% names(gl)))
  au <- augment(fit)
  expect_s3_class(au, "tbl_df")
  expect_equal(nrow(au), 96L)
  expect_true(all(c("y", "time", "arm", "id", ".fitted") %in% names(au)))
  expect_true(all(au$.fitted > 0 & au$.fitted < 1))
  nd <- tibble::tibble(time = c(0, 1), arm = c(0, 1))
  expect_equal(nrow(augment(fit, newdata = nd)), 2L)
})

test_that("fits work inside dplyr and purrr workflows", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  d <- dplyr::bind_rows(
    dplyr::mutate(fc_test_data(20, seed = 1), site = "A"),
    dplyr::mutate(fc_test_data(20, seed = 2), site = "B")
  )
  res <- d |>
    dplyr::group_by(site) |>
    dplyr::group_modify(~ tidy(fit_fewclust(.x, y ~ time + arm, id = id, time = time))) |>
    dplyr::ungroup()
  expect_equal(nrow(res), 6L)
  fits <- purrr::map(split(d, d$site), \(x) fit_fewclust(x, y ~ time + arm, id = id,
                                                         time = time))
  expect_length(fits, 2L)
  expect_true(all(purrr::map_lgl(fits, \(f) inherits(f, "fewclust"))))
})
