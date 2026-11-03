# The Vedovelli procedure for marginalized multilevel models (MMM): implementation spec

> Status: specification, to be implemented in `fewclust` **after the paper** (PI, 2026-10-21).
>
> **Note for fewclust 0.3.0 (2026-11-02).** The rule is implemented in `R/mmm.R` (`fit_mmm()`, without
> depending on MMLB), with the record asked in section 3 (`$model$check$path`, `$reason`, `$trend`,
> `$constant`) and tests on the fixtures of section 7 (`tests/testthat/test-mmm-check.R`). Section 1 below
> describes two failure modes found in S6; the final paper (manuscript v5, Section 2.4 and Supplement S16)
> describes them as **non-convergence** and **convergence with a collapsed standard error**, and no longer
> states the "shifted estimate" mode. The package documentation follows the paper. The rule itself
> (threshold 1e-3, nlm codes other than 1 and 2, fallback to `lv.formula = ~ 1`, t with G - p) is unchanged.
> Results of the 0.3.0 checks: `dev/aggiornamento_0.3.0.md`.
> Nickname chosen by the PI: *procedura Vedovelli per gli MMM*. In code and docs use a
> neutral name, e.g. `mmm_check = TRUE` / "MMM with convergence check".
> Source of truth for the rule and its evidence: branch `morris-v4`, decisions D-32, D-33,
> D-34, D-35 and step S6b (`docs/steps/S6b.md`, `reports/morris/S6b.html`).

## 1. What problem it solves

`MMLB::mm()` models the latent (between-cluster) standard deviation through `lv.formula`, and
**exponentiates every coefficient**: with `lv.formula = ~ time`,
`sigma_t = exp(a0) + exp(a1) * time`. A slope of the latent SD can therefore only be positive,
and zero is the boundary of the parameter space (`a1 -> -Inf`).

When the true slope is small (in our study about 0.01 per occasion), `nlm` often runs into that
boundary. The S6 diagnostics (`morris-v4:reports/morris/S6.html`, section 3) found two failure
modes:

1. **shifted estimate** (nlm code 5, maximum step size exceeded): the SE is right, but the
   estimate of the mean coefficients moves away from the maximum by an amount that does not
   shrink with the SE;
2. **collapsed SE** (nlm codes 2-3, slope at about 0): the estimate is right, but its SE is less
   than half of the SE of the fit without slope on the same data.

Code 4 (iteration limit, 100 by default) is harmless by itself: refitted with 2,000 iterations
those fits converge after a median of 105-140 iterations and the estimate does not change
(`morris-v4:tables/V4_S6_mmm_iterlim.csv`). The validated rule nevertheless treats code 4 as
"at the boundary" (section 2).

The type I error of the Wald test then **grows** with the number of clusters and of occasions
(in the study up to about 10% with 8-10 occasions; 138 of 864 null design points clearly above
7.5% before the fix).

## 2. The procedure (as validated)

Input: a binary outcome, cluster id, mean model `f_mean`, latent-variability model `f_lv` with at
least one non-intercept term, transition model `f_t`, quadrature points `q`.

1. **Fit A**: `mm(f_mean, lv.formula = f_lv, t.formula = f_t, data, id, q = 20)` with MMLB
   defaults. Capture warnings (do not drop them).
2. **Boundary check** — fit A is *at the boundary* if either holds:
   - `fit$control["convergence_code"]` is **not** 1 or 2 (nlm codes; 1 = gradient close to
     zero, 2 = successive iterates within tolerance; 3 = no lower point found, 4 = iteration
     limit, 5 = maximum step size exceeded five times);
   - `min_k exp(a_k) < tau`, over the latent-SD coefficients **other than the intercept**
     (names `log(sigma):<term>` in `fit$alpha`, excluding `log(sigma):(Intercept)`), with
     **`tau = 1e-3`**.
3. **Fit B** (only if A is at the boundary): the same call with **`lv.formula = ~ 1`**
   (constant latent SD), same `f_mean`, same `f_t`, same `q`. With several latent-SD terms
   (e.g. `~ time + dose`) **all** non-intercept terms are dropped together.
4. **Choice**: use fit B if it succeeded, contains the tested term, and gives a finite estimate
   and SE; otherwise keep fit A ("fallback failed").
5. **Inference**: Wald statistic `est / se` with a **t reference on `G - p` degrees of freedom**
   (`G` clusters, `p` mean-model coefficients including the intercept; `p = 4` for
   `~ time * trt`). Do **not** use the normal p-value printed by `summary.MMLB()`: with the
   normal reference the corrected MMM is clearly above 7.5% in 4 of 1,896 design points, with
   the t reference in none.
6. If `f_lv` has no non-intercept term, the procedure is a no-op (steps 2-4 are skipped).

The rule is fixed in advance: it must never be switched on or off after looking at the result.

Reference code (study engine, base R + data.table): `morris-v4:R/s6/s6_mmm_cc_run.R`
(functions `fitw`, `slope_min`, block inside `run_task`) and `morris-v4:R/s6b/s6b_ns_common.R`
(`cc_fits`). These are the exact implementations behind every number below.

## 3. What the fit object must record

Every result of the procedure must carry, inspectable and printable:

- `path`: `"A_kept"` (not at boundary), `"B_used"`, `"B_failed_A_kept"`;
- `reason` when at the boundary: `"nonconvergence (nlm code k)"` and/or
  `"latent-SD slope at boundary (min exp(a) = x < tau)"`;
- estimate, SE and p (t reference) of **both** fits when both were computed, plus warnings of each;
- `tau`, `q`, the degrees of freedom used.

Users must be able to report how many fits fell back and why (the paper will ask them to).

## 4. Evidence (all numbers from `morris-v4` tables; do not retype, re-read if needed)

| claim | value | source |
|---|---|---|
| study reanalysis: usable MMM fits in the scenarios with a latent-SD slope | 1,685,204 | `tables/V4_S6cc_gates.csv`, `tables/V4_S6b_thr_inventory.csv` |
| fit B used with `tau = 1e-3` | 1,139,699 | `tables/V4_S6b_thr_inventory.csv` |
| null design points clearly above 7.5% (of 864): study fit -> with the check | 138 -> 0 | `tables/V4_S6cc_preds.csv` (P10) |
| `tau` sensitivity (1e-4, 1e-2 vs 1e-3): mean type I / power differences | at most 0.0012 points; 0 of 156 manuscript sub-claims change | `tables/V4_S6b_thr_preds.csv` |
| new scenarios with a large true slope (latent SD 1 -> 1.9 or 3.25; random-slope variance x10, x40), 768,000 datasets: points clearly above 7.5% with the check and t | 0 of 384 | `tables/V4_S6b_ns_preds.csv` (PB2) |
| same, without the check | 8 of 384 | `tables/V4_S6b_ns_rule_safety_dp.csv` |
| share of fits where B is used, null arm, new scenarios | 35.6-51.4% (mostly nlm non-convergence with 2-3 occasions or few clusters) | `tables/V4_S6b_ns_preds.csv` (PB1), `tables/V4_S6b_ns_diag_fallback.csv` |
| type I among fits where B was used (null arm) | 5.00-5.23% | `tables/V4_S6b_ns_fallback.csv` |
| all 1,896 design points (study + new), MMM with check and t | 0 clearly above 7.5%, 77 clearly below 2.5%, answers 91.0%, power per attempt 21.2% | `tables/V4_S6b_rule_all.csv` |

Conclusion used in the paper (D-35): MMM with this procedure and the t reference is the
recommended strategy at any number of clusters and occasions; GEE with the Fay-Graubard
correction and t is the alternative.

## 5. Validated range and open questions (do not claim beyond this)

Validated: binary outcome, logit link, `MMLB::mm` defaults with `q = 20`, 10-145 clusters, 2-10
occasions, 70% complete series (30% truncated), mean model `~ time * trt` (+ one covariate in
M5), `lv.formula` `~ time` or `~ time + dose`, `t.formula` `~ 1`, `~ trt`, `~ trt + steroids`;
data generated by GenBinaryY, by a GLMM with random intercept and slope, and by a GLMM with AR(1).

Not tested (decide before implementing, or document as limits):

1. dropping only the offending latent-SD term instead of all of them;
2. the robust (sandwich) SE of `mm()` as an alternative to the fallback (in S6 the robust SE did
   **not** fix the excess; do not offer it as a substitute);
3. other `q` or other optimizers; for code 4 only, raising `iterlim` instead of falling back
   (harmless in S6, see section 1, but the validated rule falls back);
4. how to treat the H0 truth when data come from a GLMM with a large random slope (the marginal
   null is then only approximate, up to -0.008 here);
5. the behaviour of the check when `nlm` returns code 2 with a slope just above `tau`
   (harmless in our data: estimates with and without slope agree to about 1e-3 SE).

## 6. Implementation checklist (for later)

- [ ] add `MMLB` to `Suggests` (installed from GitHub `mercaldo/MMLB`; not on CRAN — check CRAN
      policy before making it a hard dependency);
- [ ] engine `fit_mmm()` with arguments `check = TRUE`, `tau = 1e-3`, `q = 20`, reference `"t"`;
- [ ] record `path`, `reason`, both fits (section 3); print and `summary` methods show them;
- [ ] decision path: add the MMM strategy as the recommended one (D-35), FG with t as alternative;
      keep the warning for 2-3 occasions and up to 25 clusters (no method is both available and
      calibrated there);
- [ ] unit tests on the fixtures in `fewclust/dev/vedovelli_fixtures/` (section 7): every branch,
      estimates to 1e-8, `path` and `reason` correct, t p-values on `G - 4` df;
- [ ] test that the procedure is a no-op when `lv.formula = ~ 1`;
- [ ] process-level time limit around `mm()` (as for the other engines): `nlm` can be slow on
      sparse data.

## 7. Fixtures

`fewclust/dev/vedovelli_fixtures/`: study datasets (scenario M3, GenBinaryY, 25 clusters, 8
occasions, null arm) regenerated with their seed, two per branch of the rule (`keep_converged`, `fallback_nonconvergence`, `fallback_slope_at_boundary`, and
the fourth branch, `fallback_failed_keep_study_fit`, does not occur in that cell; build a
synthetic test for it, e.g. by mocking a failing fit B), with `expected.csv`: stored fit A
(estimate, SE, nlm code, minimum latent-SD slope), branch, the estimate/SE used, the normal and
t p-values. The script that wrote them refits A and B and checks the stored values
(`morris-v4:R/s6b/s6b_vedovelli_fixtures.R`). Model call to reproduce them:
`mm(binY ~ time * trt, lv.formula = ~ time, t.formula = ~ 1, data, id = id, q = 20)`.

Fixture `case_04` shows the collapsed-SE mode: nlm code 3, latent-SD slope estimated at 7.3,
SE 0.099 for the interaction against 0.225 without the slope. `case_03` shows the shifted
estimate: code 5, estimate -0.090 against -0.073 without the slope, similar SE.
