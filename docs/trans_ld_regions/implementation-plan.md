# Trans LD regions implementation plan

Goal: deliver the approved Terra WDL and its reproducible runtime bundle.
Spec: design.md. Execution: inline, in an isolated prepare_QTL worktree; preserve existing analysis behavior.

- [x] Implement association and ancestry validation. Tests must reject duplicate assignments, absent chromosomes, bad coordinates and unresolved gs:// inputs; preserve threshold-passing association rows.
- [x] Implement ancestry LD tasks. Test strict MAF, allele-count diagnostics, different ancestry extents, adaptive search, missing seeds and insufficient samples. Real PLINK smoke tests must exercise BED and PGEN inputs.
- [x] Implement conservative interval merging. Hand-computed tests must cover earliest/latest union, padding, overlaps, fallback, chromosome ends, source evidence and preservation of phenotype links.
- [x] Implement WDL with explicit File inputs and command-time file lists. Run miniwdl and Cromwell womtool validation, an AST regression for workflow file writes, and rendered-command cloud-path localization tests for all three tasks.
- [x] Add micromamba Dockerfile, pinned environment, example Terra inputs, README, and GitHub Actions validation and image smoke tests. Do not publish or run Terra jobs.
- [x] Review changes and run final tests. Report local versus CI versus Terra validation separately.

Completed locally on 2026-09-07. CI build and Terra execution were not performed; see validation.md.
