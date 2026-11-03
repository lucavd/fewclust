# Every error, warning, message and print note says what happens and what to
# do, in plain words, without quoting internal counts of the simulation study.
forbidden <- "design point|Table [0-9]|scenario|2,376|mechanism|of the study|Section [0-9]"

catch_all <- function(expr) {
  out <- character(0)
  res <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      out <<- c(out, conditionMessage(e))
      NULL
    }),
    warning = function(w) {
      out <<- c(out, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      out <<- c(out, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  if (!is.null(res)) out <- c(out, utils::capture.output(print(res)),
                               utils::capture.output(print(summary(res))))
  out
}

test_that("errors are classed and say what to do", {
  d <- fc_test_data(20, seed = 1)
  const <- d
  const$y <- 0L
  expect_error(fit_fewclust(const, y ~ time + arm, id = id), "always 0",
               class = "fewclust_error")
  one <- d
  one$y <- 0L
  one$y[one$id == levels(one$id)[1L]] <- c(1L, 0L, 1L, 0L)
  expect_error(fit_fewclust(one, y ~ time + arm, id = id), "only one group",
               class = "fewclust_error")
  sep <- d
  sep$w <- rnorm(nrow(sep))
  sep$y <- as.integer(sep$w > 0)
  expect_error(fit_fewclust(sep, y ~ w + arm, id = id), "separat",
               class = "fewclust_error")
  big <- fc_test_data(10, seed = 1, n_occasions = 6)
  for (k in 1:9) big[[paste0("v", k)]] <- rnorm(nrow(big))
  expect_error(fit_fewclust(big, y ~ time + v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9,
                            id = id), "more groups than coefficients",
               class = "fewclust_error")
  single <- fc_test_data(30, seed = 2, n_occasions = 1)
  expect_error(fit_fewclust(single, y ~ arm + x, id = id), "glm\\(\\)",
               class = "fewclust_error")
  bin <- d
  bin$y <- bin$y + 1L
  expect_error(fit_fewclust(bin, y ~ arm, id = id), "binary", class = "fewclust_error")
})

test_that("warnings fire where the results are less reliable", {
  few <- fc_test_data(8, seed = 5)
  expect_warning(fit_fewclust(few, y ~ time + arm, id = id, time = time),
                 "Only 8 groups", class = "fewclust_warning")
  two <- fc_test_data(12, seed = 4, n_occasions = 2)
  expect_warning(fit_fewclust(two, y ~ time + arm, id = id, time = time, method = "gee"),
                 "2 occasions and fewer than 15 groups", class = "fewclust_warning")
  d <- fc_test_data(30, seed = 2)
  expect_warning(fit_fewclust(d, y ~ time + arm, id = id, time = time, method = "gee",
                              correction = "none"),
                 "too small", class = "fewclust_warning")
  expect_warning(fit_fewclust(d, y ~ time + arm, id = id, time = time, method = "gee",
                              correction = "kc"),
                 "fg", class = "fewclust_warning")
  holed <- d[!(d$id == "1" & d$time == sort(unique(d$time))[2L]), ]
  expect_warning(fit_fewclust(holed, y ~ time + arm, id = id, time = time),
                 "miss an occasion", class = "fewclust_warning")
})

test_that("messages are informative and can be silenced", {
  d <- fc_test_data(30, seed = 2)
  expect_no_message(fit_fewclust(d, y ~ time + arm, id = id, time = time))
  local_mocked_bindings(fit_mmm = function(...) fewclust:::fc_failure("too few events"),
                        .package = "fewclust")
  expect_message(fit_fewclust(d, y ~ time + arm, id = id, time = time),
                 "GEE", class = "fewclust_message")
})

test_that("print notes cover the small designs", {
  small <- fc_test_data(20, seed = 1, n_occasions = 3)
  fit <- fit_fewclust(small, y ~ time + arm, id = id, time = time)
  expect_match(paste(fit$notes, collapse = " "), "flip_test\\(\\)")
  # seed 1: the marginalized model gives a result (with seed 4, used up to
  # 0.4.0, its refit without the trend does not converge: no result, GEE)
  ten <- fc_test_data(10, seed = 1, n_occasions = 4)
  fit10 <- fit_fewclust(ten, y ~ time + arm, id = id, time = time)
  expect_equal(fit10$method, "mmm")
  expect_match(paste(fit10$notes, collapse = " "), "conservative")
  big <- fc_test_data(40, seed = 1, n_occasions = 6)
  expect_length(fit_fewclust(big, y ~ time + arm, id = id, time = time)$notes, 0L)
})

test_that("no message quotes the internal counts of the simulation study", {
  d <- fc_test_data(30, seed = 2)
  two <- fc_test_data(12, seed = 4, n_occasions = 2)
  few <- fc_test_data(8, seed = 5)
  small <- fc_test_data(20, seed = 1, n_occasions = 3)
  const <- d
  const$y <- 1L
  texts <- c(
    catch_all(fit_fewclust(d, y ~ time * arm, id = id, time = time)),
    catch_all(fit_fewclust(d, y ~ time * arm, id = id)),
    catch_all(fit_fewclust(two, y ~ time * arm, id = id, time = time)),
    catch_all(fit_fewclust(two, y ~ time * arm, id = id, time = time, method = "gee")),
    catch_all(fit_fewclust(few, y ~ time * arm, id = id, time = time)),
    catch_all(fit_fewclust(small, y ~ time * arm, id = id, time = time)),
    catch_all(fit_fewclust(const, y ~ time * arm, id = id, time = time)),
    catch_all(fit_fewclust(d, y ~ time * arm, id = id, time = time, method = "gee",
                           correction = "md")),
    catch_all(print(flip_test(const, y ~ time * arm, id = id, seed = 1))),
    catch_all(print(flip_test(d, y ~ time * arm, id = id, seed = 1)))
  )
  expect_gt(length(texts), 20L)
  expect_false(any(grepl(forbidden, texts)), info = paste(texts[grepl(forbidden, texts)],
                                                           collapse = "\n"))
})
