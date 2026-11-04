# fewclust 0.4.3

Preparation of the first CRAN release; the fits and their numbers are
those of 0.4.2.

- The notes printed with a fit are a clean bulleted list.
- The header of `R/constants.R` records the checksums of the final tables
  of the study; the stored values are unchanged.
- `DESCRIPTION` names the copyright holder; README shows the CRAN
  installation.

# fewclust 0.4.2

## Numbers of the simulation study

The numbers of the simulation study stored in the package (shown in
`vignette("simulation-study")`) now come from the final tables of the study
with its current rules: the refit of the convergence check of the
marginalized model counts only if it converges (the rule implemented in
0.4.1), and datasets with separation stay results.  Only the rows of the
marginalized model changed: its mean share of usable datasets is lower
(overall and with short series and few clusters), it is clearly below 2.5%
at more design points (also with 10 clusters), its unconditional power is
slightly lower, and its distances from 5% in the head-to-head map with GEE
AR(1) and Fay-Graubard change (two cells change class: 15-25 clusters with
3 occasions is now about the margin, more than 25 clusters with 2 occasions
is now similar).  It is still never clearly above 7.5%.  The values are in
the vignette; the list of changes, old and new, is in
`dev/checks_0.4.2/out/constants_diff.csv` of the source repository.  Every stored value now names, in `R/constants.R`, its source
table and the identifier of the study's registry of numbers.  The
statements of the vignette that depend on these numbers are computed from
them.  Fitting functions and their results are unchanged.

# fewclust 0.4.1

## Convergence check of the marginalized model

The convergence check now follows the rule adopted by the simulation study
behind the package (decision of 3 November 2026).  A fit with the time
trend of the latent standard deviation that is at the boundary (`nlm()`
code other than 1 or 2, a slope of the latent standard deviation below
`check_threshold`, or, with `check_coherence = TRUE`, a log-likelihood
below that of the fit without the trend) is replaced by the fit without
the trend **only if `nlm()` reports convergence of that refit (code 1 or
2) and its estimates and standard errors are finite**.  Otherwise the
marginalized model gives no result for the dataset: with
`method = "auto"` `fit_fewclust()` falls back to GEE (AR(1) working
correlation, Fay-Graubard correction, as for any other failure of the
marginalized model) with a message, and with `method = "mmm"` it stops with
an error of class `"fewclust_error"`.  The message says why: the refit did
not converge (with its `nlm()` code), or it gave no finite estimate.

Earlier behaviour (0.3.x, 0.4.0): the refit was used whenever its estimates
and standard errors were finite, whatever its `nlm()` code, and when it was
not finite the fit at the boundary was reported.  The study judged that a
fallback fit must pass the same convergence test as the fit it replaces.
Results change only for the datasets on which the check intervenes and the
refit does not converge or is not finite; on every other dataset estimates,
standard errors, statistics, degrees of freedom and p-values are identical
to 0.4.0.

- The record of the check (`fit$convergence`) has the new path
  `"no result"` and a field `refit_status` (`"converged"`,
  `"not converged"`, `"not finite"`).  The path
  `"constant fit failed, trend fit kept"` no longer exists, and
  `"constant fit used"` now always means a converged refit.  When GEE is
  reported because the marginalized model gave no result,
  `fit$convergence` keeps the record of its check and `summary()` prints
  it; the error of `method = "mmm"` carries it as `convergence`.
- The note "the reported fit's optimiser returned nlm code ..." is gone:
  with the check on, the reported fit of the marginalized model has always
  converged.
- The constants of the simulation study shown in
  `vignette("simulation-study")` are unchanged in this version.

## Agreement with the study

On 240 datasets of the study (scenario M3, the model the package fits)
where the study's rule gives no result, the package gives none on 238
(191 by the convergence check, 47 because `fit_fewclust()` refuses the
design); on the other 2 its optimiser converges where `MMLB::mm` does not,
so the check does not intervene.  On 240 datasets where the study uses a
converged refit, the package reports the same estimate and standard error
on 224; the 16 others are design refusals (7), refits with a non-finite
standard error in the package (3), trend fits that converge in the package
(3) and standard errors that differ (3).  Details in
`dev/aggiornamento_0.4.1.md` of the source repository.

# fewclust 0.4.0

fewclust is now a classic R package in the style of the tidyverse: one
fitting function with the data first, standard methods, broom verbs that
return tibbles, and short messages that say what happens and what to do.
The numbers do not change: for every method, the estimates, standard
errors, statistics, degrees of freedom and p-values are identical to those
of version 0.3.1, and so are `gee_vcov()` and `flip_test()`.

## Breaking changes

- `auto_cluster_fit(formula, data, group, time, strategy, reference, alpha)`
  is replaced by
  `fit_fewclust(data, formula, id, time, method, corstr, correction, conf_level)`.
  The correspondence:

  | 0.3.1                                   | 0.4.0                                              |
  |-----------------------------------------|----------------------------------------------------|
  | `auto_cluster_fit(f, d, group = id)`    | `fit_fewclust(d, f, id = id)`                      |
  | `strategy = "auto"`                     | `method = "auto"` (default)                        |
  | `strategy = "mmm"`                      | `method = "mmm"`                                   |
  | `strategy = "gee_ar1_fg"`               | `method = "gee"` (defaults `corstr = "ar1"`, `correction = "fg"`) |
  | `strategy = "gee_exch_md"`              | `method = "gee", corstr = "exchangeable", correction = "md"` |
  | `strategy = "gee_ind_unc"`              | `method = "gee", corstr = "independence", correction = "none"` |
  | `reference = "t"`                       | always (t on N - p degrees of freedom)             |
  | `alpha = 0.05`                          | `conf_level = 0.95`                                |
  | `fewclust_control(corstr = )`           | `corstr =` of `fit_fewclust()`                     |
  | `fit$coefficients` (data frame)         | `tidy(fit)` (tibble); `coef(fit)` is now a named vector |
  | `fit$strategy`                          | `fit$method` (`"mmm"` or `"gee"`) and `fit$fallback` |
  | `fit$model`                             | `fit$engine`                                       |
  | `fit$model$check`                       | `fit$convergence`                                  |
  | `flip_test(f, d, group = id, B = 1999)` | `flip_test(d, f, id = id, reps = 1999)`            |

- A design that cannot give an estimate is now an error of class
  `"fewclust_error"` (it was a `"fewclust_refusal"` object).  Fewer than 10
  groups is a warning, no longer a refusal.
- Removed: `cluster_strategies()`, `cluster_diagnose()`,
  `compare_strategies()`, `plot()` methods, the mixed-model strategies
  (`"glmm_wald"`, `"glmm_lrt"`), the normal reference (`reference = "z"`)
  and `fewclust_control(nAGQ = )`.  The package no longer depends on lme4,
  reformulas or graphics.
- `simulate_clustered_binary()` returns a tibble.

## New

- `id` and `time` accept bare column names or strings (tidyselect), also
  through `{{ }}` inside your own functions.
- Methods `predict()` (log odds or probability, also on `newdata`),
  `nobs()`, `formula()`, and broom's `augment()`; `tidy()` takes
  `conf.int` and `conf.level`; `glance()` gives one row per fit;
  `tidy()` also works on `flip_test()`.
- `flip_test()` warns, with the reason, when the test is not computable.

## Messages

- Errors, warnings and messages use cli and have the classes
  `"fewclust_error"`, `"fewclust_warning"` and `"fewclust_message"`.
- They no longer quote counts of the simulation study (design points,
  tables, scenarios): they say what happens and what to do.  The findings
  of the study are in `vignette("simulation-study")`.

# fewclust 0.3.1

The automatic choice now follows the corrected simulation study (manuscript
v7).  In the study a dataset whose outcome is constant (all 0 or all 1)
gives no usable result for any method; GEE with an AR(1) working
correlation and the Fay-Graubard correction had counted 9,646 such datasets
as usable.  Without them GEE with an AR(1) working
correlation, the Fay-Graubard correction and a t reference is clearly above
a 7.5% type I error at 2 design points (10 clusters, 2 occasions), and the
rule of the study selects the marginalized multilevel model with the
convergence check.  The numbers the paper took from the package do not
change: `gee_vcov()` and the engine of the marginalized model compute
exactly what version 0.3.0 computed.

## Breaking changes

- The recommended strategy is the marginalized multilevel model with the
  convergence check and a t reference on N - p degrees of freedom
  (`"mmm"`), at every number of clusters and occasions.  GEE with an AR(1)
  working correlation, the Fay-Graubard correction and a t reference
  (`"gee_ar1_fg"`) is the alternative: the plan names it in
  `$alternative`, and `auto_cluster_fit()` fits it, and records it, when
  the marginalized model returns no usable estimate.  Version 0.3.0
  recommended the GEE strategy first and fell back to the marginalized
  model.  With 2-3 occasions and at most 25 clusters `flip_test()` remains
  the complement.
- When both strategies of the rule fail, the refusal gives both reasons and
  suggests `flip_test()`.

## New cautions and notes

- `"gee_ar1_fg"` (chosen by the user, or used by the fallback) with fewer
  than 15 groups and 2 occasions carries a caution with the 2 design points
  of the study where it was clearly above 7.5% (scenario M2, 10 clusters,
  2 occasions; GLMM-native and GLMM-AR1 mechanisms), with their rejections,
  usable datasets and usable share.
- The note of 0.3.0 "with 2 occasions and 15 or more groups the alternative
  was slightly closer to 5%" is gone: in the corrected direct comparison
  (Figure 1a) that difference is about the size of the margin.  A note now
  appears where one strategy was closer to 5% than the other by more than
  the margin: with 10 groups and 3 or more occasions, GEE AR(1) with
  Fay-Graubard.
- The caution about a group missing an occasion inside its series now also
  applies to the marginalized model, whose transition term links
  consecutive observations.
- The counts quoted by the cautions are those of the corrected Table 5 of
  the study (for example 2 of 2,376 design points for `"gee_ar1_fg"`, 138
  for `"gee_ar1_kc"`, 20 for `"gee_ar1_md"`); they are generated from the
  tables of the study by `dev/make_constants.R`.

## Constant outcome

- `flip_test()` returns no p-value on a dataset whose outcome is constant,
  with the reason "constant outcome" (0.3.0 said "score numerically zero").
  The signs are drawn first, so the random stream is unchanged.
  `auto_cluster_fit()` and `cluster_strategies()` already refused such a
  dataset for every strategy (rule D2).

## Agreement with the study

- On 162 datasets regenerated from the study (scenario M3, three
  mechanisms, N in {10, 25, 65}, T in {2, 4, 8}), 127 could be compared.
  The marginalized model of the package gives the estimate and standard
  error of the study in all but 2 of the 119 datasets with a usable result.
  In one (GLMM-AR1, 10 clusters, 8 occasions) neither fit converges and the
  study used the fit without the trend, whose tested coefficient is finite;
  the package requires every mean coefficient to be finite and returns no
  usable result.  In the other (GenBinaryY, 65 clusters, 4 occasions) both
  implementations stop at a false maximum of the fit with the trend, at
  different points; with the default `check_coherence = TRUE` the package
  reports the fit without the trend there.
- `check_coherence = TRUE` (default) treats as not converged a fit with a
  trend of the latent standard deviation whose log-likelihood is below that
  of the fit without the trend.  This condition is not part of the rule of
  the study; on the 127 datasets it changed the reported fit in 4.  Set
  `fewclust_control(check_coherence = FALSE)` to apply the rule of the
  study exactly.

## Unchanged

- `fewclust_control(check_coherence = TRUE)` remains the default.  The
  coherence condition is an addition of the package to the convergence
  check of the study; `check_coherence = FALSE` gives the rule of the study
  exactly.  Since the marginalized model is now the automatic choice, this
  default concerns the recommended strategy; the decision on it is open.

# fewclust 0.3.0

The package now follows the final version of the simulation study
(manuscript v5: 2,376 design points, four data-generating mechanisms
including a purely serial one).  The numbers of the paper do not change:
the sandwich corrections of `gee_vcov()` and the engine of the
marginalized multilevel model compute exactly what version 0.2.0 computed.

## Breaking changes

- The recommended strategy is GEE with an AR(1) working correlation, the
  Fay-Graubard correction and a t reference on N - p degrees of freedom
  (`"gee_ar1_fg"`), at every number of clusters and occasions (Figure 1 of
  the study).  The marginalized multilevel model with the convergence
  check and a t reference (`"mmm"`) is the alternative: the plan names it,
  and `auto_cluster_fit()` fits it, and records it, when the GEE fit returns
  no usable estimate.  Version 0.2.0 recommended the marginalized model
  first and fell back to Fay-Graubard.
- `fewclust_control(corstr = )` defaults to `"ar1"` (was `"exchangeable"`).
- GEE strategies are reported with their working correlation, as in Table 4
  of the study: `"gee_ar1_fg"`, `"gee_exch_md"`, `"gee_ind_kc"` and so on.
  The names of 0.2.0 (`"gee_fg"`, ...) are still accepted and use the
  working correlation of `fewclust_control()`; the qualified names choose
  it directly.
- The GEE is fitted as in the study: `geepack::geeglm()` with the scale
  parameter fixed at 1 (`scale.fix = TRUE`), on rows ordered by group and,
  within a group, by `time`.  With an AR(1) working correlation the order
  defines the lags; the estimated scale changed the estimate of the
  correlation parameter.  Estimates and standard errors of the GEE
  strategies can therefore differ from those of 0.2.0.

## New

- `flip_test()`: the sign-flip score test of the study (Hemerik, Goeman and
  Finos 2020) for one coefficient, with its rules of computability.  It
  returns a p-value only; the interval obtained by inverting it is not
  offered, as in the study.  On study datasets regenerated from their seeds
  it reproduces the p-values of the paper (`dev/aggiornamento_0.3.0.md`).
- With 2-3 occasions and at most 25 groups the plan names `flip_test()` as a
  complement (`$complement`), and the caution quotes the usable shares of
  Table 4 (46% for the GEE strategy, 39% for the marginalized model, 72% for
  the sign-flip test).
- With 2 occasions and 15 or more groups a note says that the alternative
  was slightly closer to 5% there, and the GEE strategy conservative.
- A caution when, with the AR(1) working correlation, a group misses an
  occasion inside its series (`cluster_diagnose()$groups_with_gaps`), or
  the values of `time` are not equally spaced.
- The fit of the marginalized model records the outcome of the convergence
  check (`$model$check$path`, `$reason`) and both fits, with and without
  the time trend (`$trend`, `$constant`), so that how often the check
  intervened can be reported.  Tests on six study datasets reproduce the
  branch, the estimate and the p-value of the study.

## Changes in messages

- All cautions quote Table 4 of the final study (2,376 design points), per
  strategy and working correlation; the counts are generated from the
  table of the paper and checked against it.  Under a normal reference the
  cautions say that the counts refer to the t reference (the final study
  gives none under the normal reference).
- The marginalized model at 10 groups: conservative, clearly below 2.5% at
  18 of the 156 design points with 10 clusters.
- The documentation of the convergence check follows the paper: two
  failure modes, non-convergence and convergence with a collapsed standard
  error.  The coherence condition (`check_coherence`) is documented as an
  addition of the package, not part of the rule of the study.
- The marginalized model reports its model-based standard error, as in the
  study; the documentation recalls that in the clinical applications of the
  study it was smaller than the spread of the estimates.

# fewclust 0.2.0

The package now encodes the decision rule of the revised simulation study
(manuscript v4, 1,896 design points).

## Breaking changes

- The recommended strategy is the marginalized multilevel model with the
  convergence check (`"mmm"`), with a t reference on N - p degrees of
  freedom, at every number of clusters and occasions.  The crossover rule
  between the mixed-model likelihood-ratio test and the Fay-Graubard
  sandwich (45/55 clusters) and the between-within degrees of freedom of
  0.1.0 are gone.
- `auto_cluster_fit()`, `cluster_strategies()` and `cluster_diagnose()`
  take `(formula, data, group, time, ...)`; the `estimand` argument is
  removed (the estimand follows from the strategy).
- `fewclust_control()` has new arguments (`q`, `check`, `check_threshold`)
  and no longer has `time_limit`, `process_guard`, `ci_method` or
  `thresholds`.
- Removed: `prepare_manuscript()`, the `gee_custom` strategy, the
  process-level time limit, the profile confidence intervals of the mixed
  model (Wald intervals are reported).
- Refusal threshold raised from 8 to 10 groups, the smallest simulated
  number of clusters.

## Changes in results

- `gee_vcov(fit, "fg")` is now the correction of Fay and Graubard (2001)
  as published and as implemented in geesmv: the multiplier
  `(1 - min(b, Q_jj))^(-1/2)` is applied also where `Q_jj` is negative.
  Version 0.1.0 truncated `Q_jj` at 0, which inflated Fay-Graubard
  standard errors by up to about a quarter with 10 clusters.

## New

- `R/mmm.R`: the marginalized multilevel model (latent variable with a
  linear trend of its standard deviation in the occasion index, first-order
  transition term, Gauss-Hermite quadrature) with the convergence check of
  the study.  On the fixture datasets it reaches the optimum of `MMLB::mm`
  along the same path: estimates agree to 1.6e-10 (relative), standard
  errors to 2.4e-7 and log-likelihoods to 7.4e-12 (`tests/testthat/fixtures`).
- The convergence check adds one coherence condition to the study's rule
  (`fewclust_control(check_coherence = TRUE)`, the default): a trend fit
  whose log-likelihood lies below that of its own sub-model without the
  trend is a false optimum and is treated as not converged.
  `check_coherence = FALSE` gives the published two-condition rule.
- The dependence model uses the occasion index (rank of `time` among its
  distinct values), so the threshold of the check is per occasion whatever
  the unit of `time`; the mean model uses `time` as given.  A mean-model
  column whose standard deviation lies outside [0.1, 10) is rescaled by a
  power of ten for the optimiser and mapped back.
- The usability criterion of the study (finite estimate, standard error
  and p-value; |estimate| < 10; 0 < SE < 10) is applied coefficient by
  coefficient: an unusable coefficient keeps its estimate and standard
  error but gets `NA` statistic, p-value and interval, with a note; a fit
  with no usable coefficient fails, and the automatic path then falls back
  to Fay-Graubard and records it.
- Cautions quote the study's counts of design points clearly above a 7.5%
  type I error, with the reference distribution they refer to.
- `time` must be distinct within a group (error otherwise); the number of
  occasions is the longest series (the design maximum of the study); a
  caution says when the distinct values of `time` are not equally spaced
  and when incomplete rows were dropped; a single observation per group,
  and an outcome that varies in a single group, are refused.

# fewclust 0.1.0

Initial release.
