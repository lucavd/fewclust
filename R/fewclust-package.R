#' fewclust: regression for clustered binary outcomes with few clusters
#'
#' Fits regression models for a binary outcome measured repeatedly within a
#' small number of groups, with tests and intervals that stay calibrated
#' with few groups and short series.
#'
#' * [fit_fewclust()] fits the model: by default the marginalized
#'   multilevel model with a convergence check and a t reference, with GEE
#'   (AR(1) working correlation, Fay-Graubard correction) as fallback.
#' * [flip_test()] gives a p-value for one coefficient with the sign-flip
#'   score test, for the smallest designs.
#' * [gee_vcov()] computes small-sample sandwich variances for a
#'   [geepack::geeglm()] fit.
#' * [simulate_clustered_binary()] and [fc_trial] provide example data.
#'
#' Fits work with the tidyverse: data first, bare column names, and the
#' broom verbs [tidy()], [glance()] and [augment()].
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats as.formula binomial coef confint glm.fit model.matrix
#' @importFrom stats model.frame model.offset model.response nlm nobs plogis
#' @importFrom stats predict formula
#' @importFrom stats pt qlogis qnorm qt rbinom rnorm setNames vcov
#' @importFrom utils globalVariables
NULL

# `.fc_id` is the cluster-id column fc_gee_model() adds to the model data; it
# is referenced through geeglm()'s non-standard evaluation of `id`.
globalVariables(".fc_id")
