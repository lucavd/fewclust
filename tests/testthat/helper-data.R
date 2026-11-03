# Shared helpers for the test suite.

fc_test_data <- function(G, seed = 1, n_occasions = 4, ...) {
  d <- simulate_clustered_binary(n_groups = G, n_occasions = n_occasions,
                                 seed = seed, ...)
  d[order(d$id), , drop = FALSE]
}

fc_test_gee <- function(G, seed = 1, corstr = "exchangeable",
                        formula = y ~ time + arm + x, ...) {
  d <- fc_test_data(G, seed = seed, ...)
  geepack::geeglm(formula, id = id, data = d,
                  family = binomial("logit"), corstr = corstr)
}

# The design measurements behind errors, warnings and notes (internal).
diag_of <- function(formula, data, group = NULL, time = NULL) {
  fewclust:::fc_diagnose(fewclust:::fc_parse(formula, data, substitute(group),
                                             substitute(time), env = parent.frame()))
}
