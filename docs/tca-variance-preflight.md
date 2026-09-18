# TCA first-update variance screening

FitTca screens variable genes before the full TCA fit. It uses the same
first-update least-squares calculation, covariates, bounds, and normalization
as TCA 1.2.1 with `vars.mle=FALSE`, `refit_W=FALSE`, and `constrain_mu=TRUE`.
Expression stays in linear CPM space. Cell-type means are bounded by the global input minimum plus
`mu_epsilon` and maximum minus `mu_epsilon`. Covariate coefficients retain
the unrestricted solver bounds. Variance MLE remains off. This does not
constrain every donor-level expression estimate; exported values are not clipped.
Constrained TCA can omit optional p-value matrices, which cleanup accepts
while still requiring all fitted parameter matrices.

- A negative or non-finite variance coefficient excludes that gene from fitting.
  This also checks the gene's contribution to the shared noise variance.
- A nonnegative coefficient below the configured lower bound is recorded as a
  warning. The gene is retained unless another coefficient requires exclusion.
- Solver errors, non-finite mean estimates, and non-finite squared residuals
  stop the task with the gene ID and calculation phase. They do not trigger
  automatic gene removal.
- High expression alone is not an exclusion rule.

The check adds one mean/variance estimation pass before the full fit. It runs
serially and does not change the user's parallel setting for the subsequent fit.
It screens the first update only. Later iterations can still fail, and removing
genes changes the shared noise estimate in subsequent fitting. Completion of
this check does not establish convergence or validate downstream QTL results.

## Reports and model alignment

- The existing `tca_excluded_genes.tsv` lists constant genes and preflight
  removals, with reason `invalid_first_update_variance` for the latter.
- The existing `tca_model.log` records each flagged gene, component, numerical
  estimate, lower bound, reason, and action. Its summary counts removed genes
  and warning coefficients separately.
- Saved models retain the coefficient-level table in
  `model$tca_preflight$report`, plus the variable-gene order and preflight
  exclusions. This record is an audit table, not a fitted parameter matrix.
- CleanTcaModel combines preflight and post-fit exclusions in the existing
  `gene_filter` record. The unchanged export code uses that record to select
  genes from the source BED in model order. Model restart retains the record.

The fit stage includes the preflight helper. Shared fit and cleanup helpers also
require refreshed export and downstream images under the stage dependency policy.
There are no new WDL inputs or outputs. Fit-stage runtime tests cover constrained
mean bounds and cleanup with absent optional p-values. Local tests do not establish
Terra compatibility; the complete workflow has not been tested on Terra with this change.
