# Fixtures for the MMM convergence check ("Vedovelli procedure")

Six study datasets (scenario M3, GenBinaryY generator, 25 clusters, 8 occasions, null arm),
regenerated with their seed `20260805 + 1000*3 + sim_id` (L'Ecuyer-CMRG) by
`morris-v4:R/s6b/s6b_vedovelli_fixtures.R`. Columns: `id`, `time`, `trt`, `binY`.

`expected.csv` holds, per file: the stored study fit A (estimate, SE, nlm code, minimum
latent-SD slope exp(a)), whether A is at the boundary, which fit is used (A or B), the estimate
and SE used, the normal p-value printed by MMLB and the p-value with the t reference on
G - 4 = 21 degrees of freedom. `refit_A_matches` / `refit_B_matches`: refits agree with the
stored values to 1e-8.

Model: `MMLB::mm(binY ~ time * trt, lv.formula = ~ time, t.formula = ~ 1, data, id = id, q = 20)`;
fit B uses `lv.formula = ~ 1`. Spec: `fewclust/dev/vedovelli_procedure_mmm.md`.
