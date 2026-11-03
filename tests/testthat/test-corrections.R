# Empirical characterization on balanced designs: on this grid every
# correction inflates every standard error relative to the uncorrected
# sandwich.  This is NOT a per-sample theorem (see the violation test
# below); what is pinned here is (a) the aggregate behaviour on balanced
# data and (b) the absence of the garbage-values regime observed in the
# reference implementation.
test_that("on balanced designs the corrections inflate the SE (empirical)", {
  grid <- expand.grid(G = c(10, 20, 40), seed = c(1, 2, 3),
                      corstr = c("independence", "exchangeable"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    fit <- fc_test_gee(grid$G[i], seed = grid$seed[i],
                       corstr = grid$corstr[i])
    se_unc <- sqrt(diag(gee_vcov(fit, "unc")))
    for (type in c("fg", "md", "kc")) {
      se_c <- sqrt(diag(gee_vcov(fit, type)))
      expect_true(
        all(se_c >= se_unc - 1e-10),
        label = sprintf("%s >= unc at G=%d seed=%d corstr=%s",
                        type, grid$G[i], grid$seed[i], grid$corstr[i])
      )
    }
  }
})

test_that("uncorrected sandwich reproduces geepack's own robust variance", {
  fit <- fc_test_gee(30, seed = 7)
  ours <- gee_vcov(fit, "unc")
  theirs <- fit$geese$vbeta
  expect_equal(unname(ours), unname(as.matrix(theirs)), tolerance = 1e-6)
})

test_that("corrections vanish as the number of groups grows", {
  fit <- fc_test_gee(200, seed = 5)
  se_unc <- sqrt(diag(gee_vcov(fit, "unc")))
  for (type in c("fg", "md", "kc")) {
    se_c <- sqrt(diag(gee_vcov(fit, type)))
    expect_lt(max(se_c / se_unc), 1.02)
  }
})

test_that("KC and MD fail loudly, with the cluster named, when leverage hits 1", {
  # a covariate that identifies a single cluster puts that cluster's
  # leverage at 1: the corrections are undefined, and the honest answer is
  # a structured failure, not a crash or an unbounded loop
  d <- fc_test_data(6, seed = 9)
  d$w <- as.integer(d$id == levels(d$id)[1L])
  d$y[d$id == levels(d$id)[1L]] <- c(0L, 1L, 0L, 1L)  # mixed: no separation
  fit <- geepack::geeglm(y ~ time + w, id = id, data = d,
                         family = binomial("logit"), corstr = "independence")
  expect_error(gee_vcov(fit, "kc"), class = "fewclust_correction_error")
  expect_error(gee_vcov(fit, "kc"), "leverage")
  expect_error(gee_vcov(fit, "md"), class = "fewclust_correction_error")
})

test_that("gee_vcov validates its input", {
  expect_error(gee_vcov(lm(y ~ x, data = fc_trial)), "geeglm")
  fit <- fc_test_gee(10, seed = 3)
  expect_error(gee_vcov(fit, "fg", fg_bound = 1), "strictly inside")
  expect_error(fewclust_control(fg_bound = 1), "strictly inside")
})

test_that("prior weights are honoured exactly as geepack honours them", {
  d <- fc_test_data(30, seed = 3)
  d$w <- rep(c(1, 3), length.out = nrow(d))
  m <- geepack::geeglm(y ~ time + arm, id = id, data = d, weights = w,
                       family = binomial("logit"), corstr = "exchangeable")
  expect_equal(unname(gee_vcov(m, "unc")),
               unname(as.matrix(m$geese$vbeta)), tolerance = 1e-6)
})

test_that("non-contiguous ids are blocked exactly as geepack blocks them", {
  d <- fc_test_data(15, seed = 4)
  d <- d[order(d$time, d$id), ]  # geepack now sees many singleton runs
  m <- geepack::geeglm(y ~ time + arm, id = id, data = d,
                       family = binomial("logit"), corstr = "independence")
  expect_equal(unname(gee_vcov(m, "unc")),
               unname(as.matrix(m$geese$vbeta)), tolerance = 1e-6)
})

test_that("fits using geeglm's waves argument are rejected", {
  d <- fc_test_data(20, seed = 5)
  d$wave <- as.integer(round(d$time * 3)) + 1L
  m <- geepack::geeglm(y ~ time + arm, id = id, data = d, waves = wave,
                       family = binomial("logit"), corstr = "ar1")
  expect_error(gee_vcov(m, "unc"), class = "fewclust_correction_error")
  expect_error(gee_vcov(m, "unc"), "waves")
})

test_that("a garbage variance is a classed error, never a returned number", {
  # two occasions, ten groups: MD/KC leverage slips past the early guards
  # and used to come back as astronomically negative diagonals
  fit <- suppressWarnings(fc_test_gee(10, seed = 10, n_occasions = 2))
  for (type in c("md", "kc")) {
    V <- tryCatch(gee_vcov(fit, type), error = function(e) e)
    if (inherits(V, "error")) {
      expect_s3_class(V, "fewclust_correction_error")
    } else {
      expect_true(all(is.finite(V)) && all(diag(V) > 0))
    }
  }
})

test_that("Fay-Graubard follows Fay and Graubard (2001) also where diag(Q) is negative", {
  fg_reference <- function(fit, b = 0.75) {
    X <- model.matrix(fit); y <- fit$y; mu <- fit$fitted.values
    id <- as.character(fit$id); a <- fit$geese$alpha
    cl <- split(seq_along(y), factor(id, levels = unique(id)))
    parts <- lapply(cl, function(i) {
      n <- length(i); R <- matrix(a, n, n); diag(R) <- 1
      A <- sqrt(mu[i] * (1 - mu[i])); Vi <- (A %o% A) * R
      D <- X[i, , drop = FALSE] * (mu[i] * (1 - mu[i]))
      list(DtVi = crossprod(D, solve(Vi)), D = D, e = y[i] - mu[i])
    })
    B <- Reduce(`+`, lapply(parts, function(p) p$DtVi %*% p$D)); Binv <- solve(B)
    M <- Reduce(`+`, lapply(parts, function(p) {
      Q <- p$DtVi %*% p$D %*% Binv
      H <- diag((1 - pmin(b, diag(Q)))^(-0.5), nrow(Q))
      u <- H %*% p$DtVi %*% p$e; tcrossprod(u)
    }))
    list(V = Binv %*% M %*% Binv,
         any_negative = any(sapply(parts, function(p) any(diag(p$DtVi %*% p$D %*% Binv) < 0))))
  }
  set.seed(4)
  d <- fc_test_data(12, seed = 4)
  d <- d[-sample(nrow(d), 10), ]           # unbalanced clusters
  m <- geepack::geeglm(y ~ time * arm + x, id = id, data = d,
                       family = binomial("logit"), corstr = "exchangeable")
  ref <- fg_reference(m)
  expect_true(ref$any_negative)
  expect_equal(unname(gee_vcov(m, "fg")), unname(ref$V), tolerance = 1e-10)
})
