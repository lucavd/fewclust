# The engine follows the matrix notation of its references (X, Z, W, V,
# ...), so object names are not snake_case in this file.
# nolint start: object_name_linter.

#' Small-sample corrected sandwich variance for a GEE fit
#'
#' Computes the Liang-Zeger sandwich variance of a [geepack::geeglm()] fit
#' and its three small-sample corrections, implemented directly in this
#' package: Fay-Graubard (`"fg"`), Mancl-DeRouen (`"md"`) and
#' Kauermann-Carroll (`"kc"`).
#'
#' The Kauermann-Carroll factor \eqn{(I - H_i)^{-1/2}} is computed through a
#' similarity transformation that makes the matrix symmetric before the
#' eigendecomposition, so eigenvalues are real by construction.  When a
#' cluster's leverage reaches 1 (for instance when a covariate identifies a
#' single cluster) the correction is mathematically undefined: the function
#' then signals a classed error (condition class
#' `"fewclust_correction_error"`) naming the cluster, rather than returning
#' an unstable number.  A final gate rejects any non-finite or
#' non-positive-diagonal variance for every estimator.
#'
#' Prior weights and non-contiguous cluster ids are handled exactly as
#' [geepack::geeglm()] handles them (each maximal run of equal ids is one
#' cluster); fits made with the `waves` argument of [geepack::geeglm()]
#' are rejected with a classed error, since the working correlation
#' cannot be reconstructed from row positions.
#'
#' @param fit A fitted [geepack::geeglm()] object.
#' @param type `"unc"` (Liang-Zeger, uncorrected), `"fg"`, `"md"` or
#'   `"kc"`.
#' @param fg_bound Bound \eqn{b} for the Fay-Graubard multiplier
#'   \eqn{(1 - \min(b, Q_{jj}))^{-1/2}}; default 0.75, must lie strictly
#'   inside (0, 1).
#'
#' @return A \eqn{p \times p} variance-covariance matrix for the regression
#'   coefficients.
#'
#' @references
#' Fay MP, Graubard BI (2001). Small-sample adjustments for Wald-type tests
#' using sandwich estimators. *Biometrics* 57(4), 1198-1206.
#' Mancl LA, DeRouen TA (2001). A covariance estimator for GEE with improved
#' small-sample properties. *Biometrics* 57(1), 126-134.
#' Kauermann G, Carroll RJ (2001). A note on the efficiency of sandwich
#' covariance matrix estimation. *JASA* 96(456), 1387-1396.
#'
#' @examples
#' d <- simulate_clustered_binary(n_groups = 20, seed = 1)
#' m <- geepack::geeglm(y ~ time + arm,
#'   id = id, data = d,
#'   family = binomial("logit")
#' )
#' sqrt(diag(gee_vcov(m, "unc")))
#' sqrt(diag(gee_vcov(m, "fg")))
#' @export
gee_vcov <- function(fit, type = c("unc", "fg", "md", "kc"),
                     fg_bound = 0.75) {
  type <- match.arg(type)
  if (!inherits(fit, "geeglm")) {
    stop("`fit` must be a geeglm object", call. = FALSE)
  }
  stopifnot(
    "`fg_bound` must be a single number strictly inside (0, 1)" =
      fc_is_number(fg_bound) && fg_bound > 0 && fg_bound < 1
  )
  comp <- fc_gee_components(fit)
  V <- switch(type,
    unc = fc_meat_sandwich(comp, lapply(comp$clusters, `[[`, "U")),
    fg = fc_vcov_fg(comp, fg_bound),
    md = fc_vcov_md(comp),
    kc = fc_vcov_kc(comp)
  )
  fc_check_vcov(V, type)
  dimnames(V) <- list(comp$coef_names, comp$coef_names)
  V
}

# Final gate shared by all estimators: a variance matrix that is not
# finite with a strictly positive diagonal is never returned -- silently
# returning one is the numerically-garbage regime observed in the
# reference implementation.
fc_check_vcov <- function(V, type) {
  if (!all(is.finite(V)) || any(diag(V) <= 0)) {
    fc_correction_error(sprintf(
      paste0(
        "the \"%s\" estimator produced a non-finite or non-positive ",
        "variance (extreme cluster leverage or a degenerate design); ",
        "no number is returned"
      ),
      type
    ))
  }
  invisible(V)
}

fc_correction_error <- function(reason) {
  stop(errorCondition(reason,
    class = c(
      "fewclust_correction_error", "error",
      "condition"
    )
  ))
}

# Decompose a geeglm fit into per-cluster building blocks:
#   D_i = diag(dmu/deta) X_i,  V_i = A_i^{1/2} R_i A_i^{1/2},
#   B = sum_i D_i' V_i^{-1} D_i,  U_i = D_i' V_i^{-1} (y_i - mu_i).
# The dispersion phi multiplies V_i but cancels from every sandwich, so it
# is omitted throughout.  Prior weights enter through the variance
# function (validated against geepack's own vbeta to 1e-14); clusters are
# rebuilt exactly as geepack blocks them (maximal runs of equal id) and
# cross-checked against geese$clusz.
fc_gee_components <- function(fit) {
  X <- model.matrix(fit)
  y <- fit$y
  mu <- fit$fitted.values
  eta <- fit$linear.predictors
  fam <- fit$family
  id <- fit$id
  if (is.null(id)) {
    stop("the geeglm fit carries no cluster id", call. = FALSE)
  }
  if ("waves" %in% names(fit$call)) {
    fc_correction_error(paste0(
      "fits using geeglm's `waves` argument are not supported: the ",
      "working correlation cannot be reconstructed from row positions"
    ))
  }
  alpha <- unname(fit$geese$alpha)
  corstr <- fit$corstr

  runs <- rle(as.character(id))
  blocks <- split(seq_along(y), rep(seq_along(runs$lengths), runs$lengths))
  names(blocks) <- make.unique(runs$values)
  clusz <- unname(fit$geese$clusz)
  if (length(blocks) != length(clusz) || !all(lengths(blocks) == clusz)) {
    fc_correction_error(paste0(
      "cluster blocks reconstructed from the id vector do not match the ",
      "fit's geese$clusz; cannot compute a corrected variance"
    ))
  }

  w <- fit$prior.weights
  if (is.null(w)) w <- rep(1, length(y))
  a_vec <- fam$mu.eta(eta) # derivative of mu with respect to eta
  v_vec <- fam$variance(mu) / w # variance function, weighted

  clusters <- Map(function(idx, label) {
    n_i <- length(idx)
    D <- X[idx, , drop = FALSE] * a_vec[idx]
    R <- fc_working_cor(corstr, alpha, n_i)
    sdv <- sqrt(v_vec[idx])
    Vi <- (sdv %o% sdv) * R
    Vinv <- tryCatch(solve(Vi), error = function(e) NULL)
    if (is.null(Vinv)) {
      fc_correction_error(sprintf(
        paste0(
          "the working covariance is singular for cluster '%s': the ",
          "estimated working correlation is at a degenerate value ",
          "(e.g. exchangeable alpha at 1); try corstr = ",
          "\"independence\" or a conditional strategy"
        ),
        label
      ))
    }
    e <- y[idx] - mu[idx]
    list(
      D = D, Vi = Vi, e = e,
      DtVinv = crossprod(D, Vinv), U = crossprod(D, Vinv %*% e)
    )
  }, blocks, names(blocks))
  B <- Reduce(`+`, lapply(clusters, function(cl) cl$DtVinv %*% cl$D))
  Binv <- tryCatch(solve(B), error = function(e) NULL)
  if (is.null(Binv)) {
    fc_correction_error(paste0(
      "model-based information matrix is singular (aliased or degenerate ",
      "design, possibly separation)"
    ))
  }
  list(
    clusters = clusters, B = B, Binv = Binv, p = ncol(X),
    coef_names = colnames(X), cluster_labels = names(blocks)
  )
}

fc_working_cor <- function(corstr, alpha, n) {
  a <- if (length(alpha) >= 1L) alpha[1L] else 0
  switch(corstr,
    independence = diag(n),
    exchangeable = {
      R <- matrix(a, n, n)
      diag(R) <- 1
      R
    },
    ar1 = a^abs(outer(seq_len(n), seq_len(n), `-`)),
    fc_correction_error(paste0(
      "unsupported working correlation structure: ", corstr
    ))
  )
}

fc_meat_sandwich <- function(comp, U_list) {
  M <- Reduce(`+`, lapply(U_list, tcrossprod))
  comp$Binv %*% M %*% comp$Binv
}

# Fay-Graubard: the cluster score is rescaled by
# diag{(1 - min(b, Q_i,jj))^(-1/2)}, Q_i = D_i' V_i^{-1} D_i B^{-1}, exactly
# as in Fay and Graubard (2001) and in geesmv; needs no matrix square root.
fc_vcov_fg <- function(comp, fg_bound) {
  U_list <- lapply(comp$clusters, function(cl) {
    Q <- cl$DtVinv %*% cl$D %*% comp$Binv
    cl$U / sqrt(1 - pmin(fg_bound, diag(Q)))
  })
  fc_meat_sandwich(comp, U_list)
}

# Cluster leverage H_i = D_i B^{-1} D_i' V_i^{-1}, shared by MD and KC.
# Its symmetrised form S_i = L^{-1} D_i B^{-1} D_i' L^{-T} (V_i = L L')
# has real eigenvalues in [0, 1]; leverage at or above 1 makes both
# corrections undefined and is reported as a structured failure.
fc_vcov_md <- function(comp) {
  U_list <- lapply(seq_along(comp$clusters), function(k) {
    cl <- comp$clusters[[k]]
    H <- cl$D %*% comp$Binv %*% cl$DtVinv
    IH <- diag(nrow(H)) - H
    rc <- tryCatch(rcond(IH), error = function(err) 0)
    e_adj <- if (rc < 1e-10) {
      NULL
    } else {
      tryCatch(solve(IH, cl$e), error = function(err) NULL)
    }
    if (is.null(e_adj)) {
      fc_correction_error(sprintf(
        paste0(
          "Mancl-DeRouen correction undefined: leverage at or above ",
          "1 in cluster '%s' (I - H not invertible)"
        ),
        comp$cluster_labels[k]
      ))
    }
    cl$DtVinv %*% e_adj
  })
  fc_meat_sandwich(comp, U_list)
}

fc_vcov_kc <- function(comp) {
  U_list <- lapply(seq_along(comp$clusters), function(k) {
    cl <- comp$clusters[[k]]
    L <- tryCatch(t(chol(cl$Vi)), error = function(e) NULL)
    if (is.null(L)) {
      fc_correction_error(sprintf(
        paste0(
          "Kauermann-Carroll correction undefined: the working ",
          "covariance of cluster '%s' is not positive definite"
        ),
        comp$cluster_labels[k]
      ))
    }
    G <- solve(L, cl$D %*% comp$Binv %*% t(cl$D)) # L^{-1} D B^{-1} D'
    S <- t(solve(L, t(G))) # ... L^{-T}, symmetric
    S <- (S + t(S)) / 2
    es <- eigen(diag(nrow(S)) - S, symmetric = TRUE)
    if (min(es$values) <= 1e-8) {
      fc_correction_error(sprintf(
        paste0(
          "Kauermann-Carroll correction undefined: leverage at or ",
          "above 1 in cluster '%s' (min eigenvalue of I - S is %.2e)"
        ),
        comp$cluster_labels[k], min(es$values)
      ))
    }
    # (I - H)^{-1/2} = L Q diag(lambda^{-1/2}) Q' L^{-1}
    inner <- es$vectors %*% (t(es$vectors) / sqrt(es$values))
    e_adj <- L %*% (inner %*% solve(L, cl$e))
    cl$DtVinv %*% e_adj
  })
  fc_meat_sandwich(comp, U_list)
}
# nolint end
