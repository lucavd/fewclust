# Generates data/fc_trial.rda from the in-package generator, with a fixed
# seed, so that every example and vignette runs without external data.
# Run from the package root:
#   dev/r.sh 'source("data-raw/fc_trial.R")'
source(file.path("R", "simulate.R"))

fc_trial <- simulate_clustered_binary(
  n_groups = 24, n_occasions = 4,
  effect = 0.6, sigma_u = 1, slope_sd = 0, prevalence = 0.3,
  seed = 2026
)

if (!dir.exists("data")) dir.create("data")
save(fc_trial, file = file.path("data", "fc_trial.rda"), compress = "bzip2")
