#' Synthetic clustered binary trial
#'
#' A small synthetic dataset in the shape the package targets: a binary
#' outcome measured at 4 occasions in each of 24 groups, with a group-level
#' treatment indicator.  Generated entirely in-package by
#' `data-raw/fc_trial.R`, which calls
#' `simulate_clustered_binary(n_groups = 24, n_occasions = 4, effect = 0.6,
#' sigma_u = 1, slope_sd = 0, prevalence = 0.3, seed = 2026)`.
#'
#' @format A data frame with 96 rows and 5 variables:
#' \describe{
#'   \item{id}{group identifier (factor, 24 levels)}
#'   \item{time}{measurement occasion, scaled to \eqn{[0, 1]}}
#'   \item{arm}{group-level treatment indicator (0/1, constant within group)}
#'   \item{x}{observation-level continuous covariate}
#'   \item{y}{binary outcome (0/1)}
#' }
#' @source Simulated; see `data-raw/fc_trial.R` in the package sources.
#' @examples
#' head(fc_trial)
#' tapply(fc_trial$y, fc_trial$arm, mean)
"fc_trial"
