# Validation record — 2026-09-07

## Passed locally

- miniwdl 1.14.2: complete WDL type/syntax validation.
- Cromwell womtool 87: complete WDL validation, reported `Success!`.
- 12 unittest cases, with no skips in the final run. Used the pinned PLINK 2.0.0-a.6.9 macOS binary (29 Jan 2025), Python 3.11.4, PyArrow 19.0.1 and zstandard 0.23.0.
- Additional earlier BED/PGEN tests passed with native PLINK 2.0.0-a.6.35.
- Python and YAML syntax checks.
- Example Terra JSON checked against the workflow's typed inputs.
- Independent code review. Fixed numeric BIM chromosome X handling and the status for seeds without reported LD support; regression tests pass.

Test coverage includes:

1. BED and PGEN inputs; two ancestry groups with different LD extents.
2. Strict MAF exclusion at exactly 1%, and retention of MAF 1.5% with MAC 3.
3. Adaptive search expansion and unresolved maximum-search flags.
4. Missing variants, insufficient sample counts and constant dosage.
5. Compressed PVAR input and rejection of association/genotype position disagreement.
6. Chromosome X encoded as 23 in BIM (female diploid fixture).
7. Earliest/latest ancestry union, padding, fallback and preservation of phenotype links.
8. Minimum-width expansion at the start, middle and end of chromosomes; preservation of larger intervals.
9. Empty result output tables, Parquet/TSV source contexts and duplicate ancestry assignment rejection.
10. All three actual WDL task commands rendered after simulated cloud-to-local File localization, using real PLINK. Paths contained apostrophes, command-substitution text and backticks. Unresolved gs:// inputs failed before computation, and no injected commands ran.
11. AST regression rejecting workflow-scope file-writing functions and checking required File inputs.

Final command:

```bash
PLINK2=/path/to/plink2 python -m unittest discover -s tests/trans_ld_regions -v
```

Result: `Ran 12 tests ... OK`.

## Not run

- No local Docker image build.
- The GitHub Actions Linux image build and smoke jobs have not run for this new bundle.
- No image publication.
- No complete Terra workflow execution or real-cohort analysis.
- Mixed-sex chromosome X data have not been tested.

These local results do not constitute Terra validation. The supplied Actions jobs test the Linux runtime after the bundle is placed at a GitHub repository root. A published image and real Terra inputs are required before launch.

Repository integration: added to prepare_QTL under workflows/genotype, tools/trans_ld_regions and tests/trans_ld_regions. The additional regression checks that a newly generated cloud File remains File-typed until localization. All five repository CI path-selection tests and the complete 36-WDL file-scope check passed.
