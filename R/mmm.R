# The engine follows the matrix notation of its references (X, Z, W, V,
# ...), so object names are not snake_case in this file.
# nolint start: object_name_linter.

# Marginalized multilevel model for binary data (Heagerty 1999; Schildcrout
# and Heagerty 2007): with the convergence check and a t reference, the
# recommended strategy of the simulation study.
#
# Mean model (the estimand):   logit mu_ij = x_ij' beta
# Dependence model:            logit P(Y_ij = 1 | Y_i,j-1, b_i) =
#                                Delta_ij + gamma_ij Y_i,j-1 + sigma_ij b_i,
#                              b_i ~ N(0, 1), Y_i0 = 0,
#   with sigma_ij = Z_ij' exp(a) (so the latent standard deviation is
#   exp(a_0) + exp(a_1) time when Z = [1, time]) and gamma_ij = W_ij' g.
# Delta_ij is the intercept that makes the conditional model average back to
# mu_ij; it is solved by Newton's method, recursively along the series, as
# in MMLB::mm.  The likelihood integrates b_i by Gauss-Hermite quadrature.
#
# The latent standard deviation is written exp(a_0) + exp(a_1) time, so a
# slope of zero lies on the boundary of the parameter space.  There nlm
# fails in two ways (Section 2.4 and Supplement S16 of the study): it stops
# without converging, or it converges with a collapsed standard error; both
# produced excess rejections with long series.  The convergence check of
# the study (fit_mmm) declares the fit at the boundary when the optimiser
# does not converge or an estimated slope of the latent standard deviation
# is below `threshold`, and then replaces it by the fit with a constant
# latent standard deviation.

# Gauss-Hermite nodes and weights for a standard normal integrand
# (Golub-Welsch on the Jacobi matrix of the probabilists' Hermite
# polynomials).
fc_gh <- function(q) {
  J <- matrix(0, q, q)
  off <- sqrt(seq_len(q - 1))
  J[cbind(seq_len(q - 1), seq_len(q - 1) + 1L)] <- off
  J[cbind(seq_len(q - 1) + 1L, seq_len(q - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1L, ]^2)
}

# Arrange the design once: rows must be ordered by cluster and, within a
# cluster, by position in the series (fit_mmm orders them by `time`).  A
# column of X whose standard deviation lies outside [0.1, 10) is divided by
# a power of ten that brings it inside (nlm stops at a false optimum when a
# covariate is on a large scale); columns on the scale of the study are left
# as they are, so the optimiser follows the path of MMLB::mm.  Estimates and
# covariance are mapped back to the original scale in fc_mmm_fit_once.
fc_mmm_design <- function(y, X, id, Z, W = NULL, offset = NULL, q = 20L) {
  cl <- as.integer(id)
  pos <- stats::ave(seq_along(cl), cl, FUN = seq_along)
  ylag <- c(0L, y[-length(y)]) * (pos > 1L)
  k <- floor(log10(apply(X, 2L, stats::sd)))
  scale <- ifelse(is.finite(k) & (k < -1 | k >= 1), 10^k, 1)
  list(
    y = y, X = X / rep(scale, each = nrow(X)), scale = scale, Z = Z, W = W,
    offset = if (is.null(offset)) numeric(length(y)) else offset,
    cl = cl, n_cl = nlevels(id), ylag = ylag,
    by_pos = split(seq_along(cl), pos), gh = fc_gh(q),
    p = ncol(X), m_g = if (is.null(W)) 0L else ncol(W), m_s = ncol(Z)
  )
}

# Negative log-likelihood with its analytic gradient as attribute, in the
# form nlm() consumes.  theta = (beta, g, a).
fc_mmm_negloglik <- function(theta, d) {
  z <- d$gh$z
  w <- d$gh$w
  Q <- length(z)
  n <- length(d$y)
  npar <- length(theta)
  bad_nll <- structure(Inf, gradient = rep(NA_real_, npar))
  if (!all(is.finite(theta))) {
    return(bad_nll)
  }
  beta <- theta[seq_len(d$p)]
  g <- theta[d$p + seq_len(d$m_g)]
  a <- theta[d$p + d$m_g + seq_len(d$m_s)]
  if (any(a > 50)) {
    return(bad_nll)
  }

  eta <- drop(d$X %*% beta) + d$offset
  mu <- plogis(eta)
  gam <- if (d$m_g) drop(d$W %*% g) else numeric(n)
  sig <- drop(d$Z %*% exp(a))
  if (!all(is.finite(mu)) || !all(is.finite(gam)) || !all(is.finite(sig))) {
    return(bad_nll)
  }
  # derivatives of the three linear predictors with respect to theta
  dmu <- cbind(mu * (1 - mu) * d$X, matrix(0, n, d$m_g + d$m_s))
  dgam <- cbind(matrix(0, n, d$p), d$W, matrix(0, n, d$m_s))
  dsig <- cbind(matrix(0, n, d$p + d$m_g), d$Z * rep(exp(a), each = n))

  Delta <- numeric(n)
  dDelta <- matrix(0, n, npar)
  P <- matrix(0, d$n_cl, Q) # P(Y_i,j-1 = 1 | b = z_k)
  dP <- replicate(npar, P, simplify = FALSE)
  for (idx in d$by_pos) {
    ci <- d$cl[idx]
    Pk <- P[ci, , drop = FALSE]
    sz <- outer(sig[idx], z)
    dl <- eta[idx]
    for (iter in seq_len(100L)) {
      h0 <- plogis(dl + sz)
      h1 <- plogis(dl + gam[idx] + sz)
      A0 <- h0 * (1 - h0) * (1 - Pk)
      A1 <- h1 * (1 - h1) * Pk
      denom <- drop((A0 + A1) %*% w)
      if (!all(is.finite(denom)) || any(denom <= 0)) {
        return(bad_nll)
      }
      step <- drop(((h0 * (1 - Pk) + h1 * Pk) %*% w) - mu[idx]) / denom
      step <- pmin(0.5, pmax(-0.5, step))
      dl <- dl - step
      if (max(abs(step)) < 1e-12) break
    }
    # Recompute conditional probabilities and weights at the final dl so
    # that dDelta and dP are exact at the converged Delta.
    h0 <- plogis(dl + sz)
    h1 <- plogis(dl + gam[idx] + sz)
    A0 <- h0 * (1 - h0) * (1 - Pk)
    A1 <- h1 * (1 - h1) * Pk
    Delta[idx] <- dl
    dF_dDelta <- drop((A0 + A1) %*% w)
    if (!all(is.finite(dF_dDelta)) || any(dF_dDelta <= 0)) {
      return(bad_nll)
    }
    Az <- drop(((A0 + A1) * rep(z, each = length(idx))) %*% w)
    A1w <- drop(A1 %*% w)
    H <- h1 - h0
    for (r in seq_len(npar)) {
      dPk <- dP[[r]][ci, , drop = FALSE]
      dF <- Az * dsig[idx, r] + A1w * dgam[idx, r] +
        drop((H * dPk) %*% w) - dmu[idx, r]
      dDelta[idx, r] <- -dF / dF_dDelta
      dP[[r]][ci, ] <- A0 * (dDelta[idx, r] + outer(dsig[idx, r], z)) +
        A1 * (dDelta[idx, r] + dgam[idx, r] + outer(dsig[idx, r], z)) +
        H * dPk
    }
    P[ci, ] <- h0 * (1 - Pk) + h1 * Pk
  }

  lp <- Delta + gam * d$ylag + outer(sig, z)
  ll <- plogis((2 * d$y - 1) * lp, log.p = TRUE)
  S <- rowsum(ll, d$cl) # cluster log-likelihood at each node
  mx <- apply(S, 1L, max)
  E <- exp(S - mx) * rep(w, each = d$n_cl)
  L <- rowSums(E)
  if (!all(is.finite(L)) || any(L <= 0)) {
    return(bad_nll)
  }
  omega <- E / L # posterior weight of each node
  resid <- d$y - plogis(lp)
  grad <- vapply(seq_len(npar), function(r) {
    dlp <- dDelta[, r] + dgam[, r] * d$ylag + outer(dsig[, r], z)
    sum(omega * rowsum(resid * dlp, d$cl))
  }, numeric(1))
  nll <- -sum(log(L) + mx)
  if (!is.finite(nll) || !all(is.finite(grad))) {
    return(bad_nll)
  }
  structure(nll, gradient = -grad)
}

# Observed information by central differences of the analytic gradient.
fc_mmm_hessian <- function(theta, d, h = 1e-5) {
  g <- function(th) attr(fc_mmm_negloglik(th, d), "gradient")
  H <- vapply(seq_along(theta), function(j) {
    e <- replace(numeric(length(theta)), j, h)
    (g(theta + e) - g(theta - e)) / (2 * h)
  }, numeric(length(theta)))
  (H + t(H)) / 2
}

# Invert the observed information using Cholesky so indefinite Hessians
# (saddle points / boundary artefacts) never produce collapsed standard
# errors.  When the constant-SD model sits on the zero-variance boundary
# (sigma_0 -> 0), the log-likelihood is flat in a_0 and orthogonal to
# (beta, gamma) to O(sigma_0^2), so the non-boundary block is inverted.
fc_mmm_vcov <- function(theta, d) {
  npar <- length(theta)
  na_mat <- matrix(NA_real_, npar, npar)
  H <- tryCatch(fc_mmm_hessian(theta, d), error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) {
    return(na_mat)
  }
  V <- tryCatch(chol2inv(chol(H)), error = function(e) NULL)
  if (!is.null(V)) {
    return(V)
  }
  if (d$m_s == 1L && is.finite(theta[npar]) && exp(theta[npar]) < 1e-3) {
    sub_idx <- seq_len(npar - 1L)
    V_sub <- tryCatch(
      chol2inv(chol(H[sub_idx, sub_idx, drop = FALSE])),
      error = function(e) NULL
    )
    if (!is.null(V_sub)) {
      V <- na_mat
      V[sub_idx, sub_idx] <- V_sub
      return(V)
    }
  }
  na_mat
}

# One maximum-likelihood fit.  Starting values and optimiser settings are
# those of MMLB::mm, so that the convergence check behaves as in the study.
fc_mmm_fit_once <- function(d, inits) {
  opt <- nlm(fc_mmm_negloglik, inits,
    d = d, stepmax = 1, steptol = 1e-6,
    iterlim = 100L, check.analyticals = FALSE
  )
  theta <- opt$estimate
  V <- fc_mmm_vcov(theta, d)
  back <- c(1 / d$scale, rep(1, d$m_g + d$m_s)) # to the original scale
  theta <- theta * back
  V <- V * tcrossprod(back)
  nm <- c(
    colnames(d$X),
    if (d$m_g) paste0("gamma:", colnames(d$W)),
    paste0("log(sigma):", colnames(d$Z))
  )
  names(theta) <- nm
  dimnames(V) <- list(nm, nm)
  list(
    theta = theta, vcov = V, logLik = -opt$minimum, code = opt$code,
    iterations = opt$iterations, p = d$p,
    slope_min = if (d$m_s > 1L) {
      min(exp(theta[d$p + d$m_g + 2:d$m_s]))
    } else {
      NA_real_
    }
  )
}

fc_mmm_inits <- function(d) {
  b <- suppressWarnings(
    glm.fit(d$X, d$y, offset = d$offset, family = binomial())$coefficients
  )
  a_init <- rep(1, d$m_s)
  if (d$m_s > 1L) {
    z_max <- apply(abs(d$Z[, -1L, drop = FALSE]), 2L, max)
    # Keep MMLB::mm's starting value a = 1 on short/moderate series (<= 12);
    # scale down on longer series so exp(a) * Z never overflows plogis.
    long <- is.finite(z_max) & z_max > 12
    a_init[-1L][long] <- log(2 / z_max[long])
  }
  c(b, rep(1, d$m_g), a_init)
}

# Occasion index for the latent-SD time trend: exact grid rank on a common
# equally spaced grid (matching MMLB::mm on the study's designs), and
# elapsed time divided by the median within-group spacing on irregular or
# continuous time grids (so pooled ranks never blow up with N * T).
fc_occasion_index <- function(tv, id) {
  grid <- sort(unique(tv))
  regular <- length(grid) < 3L ||
    isTRUE(all.equal(diff(grid), rep(diff(grid)[1L], length(grid) - 1L)))
  if (regular) {
    return(match(tv, grid) - 1)
  }
  steps <- unlist(tapply(tv, id, diff), use.names = FALSE)
  step_unit <- stats::median(steps[steps > 0])
  if (!is.finite(step_unit) || step_unit <= 0) {
    return(stats::ave(seq_along(tv), id, FUN = seq_along) - 1)
  }
  (tv - min(tv)) / step_unit
}

# Status of the refit without the time trend, for the convergence check:
# "converged" (nlm code 1 or 2, finite estimates and standard errors),
# "not finite" (an error, or a non-finite estimate or standard error) or
# "not converged" (finite, nlm code 3, 4 or 5).
fc_refit_status <- function(refit, finite) {
  if (!finite) {
    return("not finite")
  }
  if (refit$code %in% c(1L, 2L)) "converged" else "not converged"
}

# Why a refit cannot replace a boundary fit, for messages.
fc_refit_problem <- function(refit, status) {
  if (inherits(refit, "fc_failure")) {
    return(paste0(
      "returned no finite estimate (",
      sub("^refit without the time trend: ", "", refit$reason), ")"
    ))
  }
  if (status == "not finite") {
    return(sprintf(
      "returned no finite estimate or standard error (nlm code %d)",
      refit$code
    ))
  }
  sprintf("did not converge (nlm code %d)", refit$code)
}

# The engine used by fit_fewclust().  With `time` the latent standard
# deviation has a linear trend in the occasion index (0, 1, 2, ... on a
# regular grid, or elapsed time scaled by the median within-group spacing
# on irregular/continuous grids, so that the threshold of the convergence
# check is per occasion whatever the unit of `time`) and the series carries
# a first-order transition term; without `time` the latent standard
# deviation is constant and there is no transition term.  The mean model
# uses `time` as given.
#
# Convergence check (the rule of the study, plus the optional coherence
# condition of fewclust_control): the fit with the trend is at the boundary
# when it fails numerically, when nlm does not report convergence (code
# other than 1 or 2), when its estimates or standard errors are not finite,
# when the estimated slope of the latent standard deviation is below
# `check_threshold`, or (coherence condition, not part of the study's rule)
# when its log-likelihood is below that of the fit without the trend, which
# is its own sub-model (a false optimum).  A boundary fit is replaced by
# the fit without the trend only when nlm reports convergence of that fit
# (code 1 or 2) and its estimates and standard errors are finite (the rule
# of the study, decision D-90); otherwise the model gives no result: an
# fc_failure carrying the record of the check, so that fit_fewclust() falls
# back to GEE (method = "auto") or stops (method = "mmm").  The boundary
# fit is never reported.  Without `time`, a non-converged or non-finite fit
# is likewise rejected when `control$check` is TRUE.
fit_mmm <- function(parsed, control, reference, alpha, df) {
  time <- parsed$time
  ord <- if (is.null(time)) {
    order(parsed$id)
  } else {
    order(parsed$id, parsed$data[[time]])
  }
  y <- parsed$y[ord]
  X <- parsed$X[ord, , drop = FALSE]
  id <- parsed$id[ord]
  offset <- parsed$offset[ord]
  one <- matrix(1, length(y), 1L, dimnames = list(NULL, "(Intercept)"))
  Z <- one
  W <- NULL
  if (!is.null(time)) {
    tv <- parsed$data[[time]][ord]
    Z <- cbind(one, occasion = fc_occasion_index(tv, id))
    W <- one
  }
  fit_once <- function(Z, what) {
    d <- fc_mmm_design(y, X, id, Z, W, offset, control$q)
    fc_quietly(fc_mmm_fit_once(d, fc_mmm_inits(d)), what)
  }
  finite <- function(f) {
    all(
      is.finite(f$theta[seq_len(f$p)]),
      is.finite(fc_se(f$vcov)[seq_len(f$p)])
    )
  }
  record <- function(f) {
    k <- seq_len(f$p)
    list(
      estimate = setNames(f$theta[k], colnames(X)),
      std.error = setNames(fc_se(f$vcov)[k], colnames(X)),
      code = f$code, slope_min = f$slope_min, logLik = f$logLik
    )
  }

  trend <- fit_once(Z, "marginalized model")
  fit <- trend$value
  trend_failed <- inherits(fit, "fc_failure")
  if (trend_failed && !(control$check && !is.null(time))) {
    return(fit)
  }
  notes <- trend$notes
  check <- list(
    requested = control$check && !is.null(time), boundary = FALSE,
    latent_sd = if (is.null(time)) "constant" else "time trend",
    path = if (is.null(time)) "no trend to check" else "trend fit kept",
    reason = character(0), threshold = control$check_threshold,
    coherence = control$check_coherence, refit_status = NA_character_,
    trend = if (is.null(time) || trend_failed) NULL else record(fit),
    constant = NULL
  )
  if (!is.null(time) && !control$check) {
    notes <- c(notes, paste0(
      "the convergence check was switched off; the ",
      "study recommends the marginalized model ",
      "only with it"
    ))
  }
  if (is.null(time) && control$check) {
    check$constant <- record(fit)
    if (!fit$code %in% c(1L, 2L)) {
      check$path <- "no result"
      check$latent_sd <- "no result"
      out <- fc_failure(sprintf(
        "the optimiser did not converge (nlm code %d)", fit$code
      ))
      out$check <- check
      return(out)
    }
    if (!finite(fit)) {
      check$path <- "no result"
      check$latent_sd <- "no result"
      out <- fc_failure(sprintf(
        "no finite estimate or standard error (nlm code %d)", fit$code
      ))
      out$check <- check
      return(out)
    }
  }
  if (check$requested) {
    constant <- fit_once(one, "refit without the time trend")
    refit <- constant$value
    ok_refit <- !inherits(refit, "fc_failure") && finite(refit)
    below <- !trend_failed && ok_refit && is.finite(fit$logLik) &&
      fit$logLik < refit$logLik - 1e-4
    why <- if (trend_failed) {
      sprintf(
        "trend fit failed (%s)",
        sub("^marginalized model: ", "", fit$reason)
      )
    } else {
      c(
        if (!fit$code %in% c(1L, 2L)) {
          sprintf("the optimiser did not converge (nlm code %d)", fit$code)
        },
        if (!finite(fit)) {
          "estimates or standard errors are not finite"
        },
        if (is.na(fit$slope_min)) {
          "the slope of the latent standard deviation is not finite"
        } else if (fit$slope_min < control$check_threshold) {
          sprintf(
            paste0(
              "the slope of the latent standard deviation is %.2g per ",
              "occasion (< %g)"
            ),
            fit$slope_min, control$check_threshold
          )
        },
        if (control$check_coherence && below) {
          paste0(
            "its log-likelihood is below that of the fit without the trend ",
            "(a false optimum)"
          )
        }
      )
    }
    check$boundary <- length(why) > 0L
    check$reason <- why
    if (!inherits(refit, "fc_failure")) check$constant <- record(refit)
    check$refit_status <- fc_refit_status(refit, ok_refit)
    if (check$boundary && check$refit_status == "converged") {
      fit <- refit
      check$latent_sd <- "constant"
      check$path <- "constant fit used"
      notes <- c(
        notes, constant$notes,
        paste0(
          "convergence check: ", paste(why, collapse = "; "),
          "; the fit without the time trend of the latent ",
          "standard deviation is reported"
        )
      )
    } else if (check$boundary) {
      check$latent_sd <- "no result"
      check$path <- "no result"
      out <- fc_failure(paste0(
        "convergence check: ", paste(why, collapse = "; "),
        "; the refit without the time trend of the latent standard ",
        "deviation ", fc_refit_problem(refit, check$refit_status),
        ", so the model gives no result"
      ))
      out$check <- check
      return(out)
    }
  }
  fit$check <- check
  k <- seq_len(fit$p)
  fc_wald_result(
    colnames(X), fit$theta[k], fc_se(fit$vcov)[k],
    fit$vcov[k, k, drop = FALSE],
    structure(fit, class = "fewclust_mmm"), reference, df, alpha,
    notes,
    scale = fc_col_scale(X)
  )
}
# nolint end
