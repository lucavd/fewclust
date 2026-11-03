# fewclust

Regression for a binary outcome measured repeatedly within a small number
of groups (clusters): 10 to about 150 groups, short series.  With few
groups the usual analyses give tests that reject too often; fewclust fits
the methods that stayed calibrated in a simulation study of these designs,
with defaults you do not need to change.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lucavd/fewclust")
```

## Usage

```r
library(fewclust)

fit <- fit_fewclust(fc_trial, y ~ time * arm, id = id, time = time)
fit
tidy(fit, conf.int = TRUE)
```

`fit_fewclust()` takes the data first, so it works in a pipe, and the
columns as bare names:

```r
fc_trial |>
  fit_fewclust(y ~ time * arm + x, id = id, time = time) |>
  tidy()
```

By default it fits the marginalized multilevel model (Heagerty 1999) with a
convergence check, with tests on a t distribution with N - p degrees of
freedom (N groups, p coefficients).  If that model gives no usable estimate,
it fits GEE with an AR(1) working correlation and the Fay-Graubard
small-sample correction instead, and tells you.  `method = "gee"` asks for
GEE directly; `corstr` and `correction` choose its variant.

The coefficients are marginal (population-averaged) log odds ratios.

## What the messages mean

- An **error** means no estimate can be computed (a constant outcome,
  separation, more coefficients than groups, ...).  It says why and what to
  do instead.
- A **warning** means the result is less reliable than usual (fewer than
  10 groups, GEE with 2 occasions and fewer than 15 groups, a correction
  other than Fay-Graubard, ...).
- A **message** reports a choice made for you (the fallback to GEE, no
  `time` given).
- **Notes**, printed with the fit, describe small designs: with 2-3
  occasions and at most 25 groups the tests have little power, and
  `flip_test()` gives a p-value that keeps its level.

All conditions have a class (`"fewclust_error"`, `"fewclust_warning"`,
`"fewclust_message"`), so they can be caught or silenced selectively.

## Functions

- `fit_fewclust()`: fit the model; methods `print()`, `summary()`,
  `coef()`, `vcov()`, `confint()`, `predict()`, `nobs()`, and the broom
  verbs `tidy()`, `glance()`, `augment()`.
- `flip_test()`: sign-flip score test of one coefficient (Hemerik, Goeman
  and Finos 2020), a p-value only, for the smallest designs.
- `gee_vcov()`: Fay-Graubard, Mancl-DeRouen and Kauermann-Carroll sandwich
  variances for a `geepack::geeglm()` fit.
- `simulate_clustered_binary()` and `fc_trial`: example data.
- `fewclust_control()`: tuning constants (rarely needed).

See `vignette("fewclust")` to get started and
`vignette("simulation-study")` for what the study found.
