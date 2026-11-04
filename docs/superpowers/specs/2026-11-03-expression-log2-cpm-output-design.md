# Expression log2-CPM output

The user approved an additional expression BED with the same genes and samples
as `.scaled`, paired with the existing scaled phenotype PCs.

## Design

Write `<OutputPrefix>.expression.log2CPM.bed.gz` from the expression preparation
task and expose it as workflow `File Log2CpmBedFile`. For count and nonnegative
CPM inputs, retain `log2(CPM + 1)` values without centering, scaling, or
residualization. Preserve an existing `Log2CpmBed` input without another log
transform. For `SignedCpm`, preserve the existing signed-log transform; those
values do not have the ordinary CPM fold-change interpretation.

Reuse the sample IDs retained by the scaled connectivity check. Use the same
annotation join and gene ordering for both BEDs. Do not run connectivity or PCA
again. Keep existing outputs and PC calculations unchanged.

## Approaches

The selected approach exports the log2 BED during preparation and pairs it with
`ScaledPhenotypePCsAllOut`. A separate log2 PCA would change gene weighting and
add a calculation. A downstream conversion would repeat work and require a
separate sample-alignment step. The user selected the additional preparation
output and existing scaled PCs.

## Outlier inputs

Use `Log2CpmBedFile` as expression `phenotype_bed` and
`ScaledPhenotypePCsAllOut` as `principal_components_tsv`. Supply fixed covariates
separately from phenotype PCs. Supply the per-ome LoF carrier table, the shared
gene annotation GTF, and the exact ome name `expression`. The default haplo drop
is `1.0`. Do not pass merged QTL covariates containing phenotype PCs or a
residualized expression BED.

## Constraints and validation

Use WDL 1.0 with typed File outputs and command logging. Keep file creation
inside the existing task. Add regression checks for input scale, retained sample
identity, gene metadata, and output discovery. Use the existing GitHub Actions
expression image tests; do not build Docker locally or submit a cloud job.
Local syntax and task tests do not establish complete Terra validation.
