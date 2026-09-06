# TCA first-update variance screening

FitTca screens variable genes before the full TCA fit. It uses the same
first-update least-squares calculation, covariates, bounds, and normalization
as TCA 1.2.1 with `vars.mle=FALSE`, `refit_W=FALSE`, and `constrain_mu=FALSE`.
Expression stays in linear CPM space. The estimator is not changed and negative
estimates are not clamped.

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

Both implementation entrypoints and the new helper are in the fit script
directory. Only the fit-stage image defaults need updating; model cleanup also
uses that stage. There are no new WDL inputs or outputs and no new GitHub checks.
Local checks cover flag classification, the reported numerical failure, solver
errors, cleanup/export alignment, and existing regression tests. A complete
Terra run and a 100-iteration fit have not been validated by these checks.
