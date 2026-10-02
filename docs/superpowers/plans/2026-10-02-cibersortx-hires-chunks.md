# CIBERSORTx HiRes gene chunks implementation plan

> For agentic workers: use the native execution path in this session, with independent boundary tests and a final review.

**Goal:** Split a supplied target-gene list into chunks of a requested maximum number of genes, run HiRes for each chunk with the full expression background and cohort, then return one merged matrix per cell type.

**Architecture:** Keep the existing `RunHiRes` task and its native parameters. Add `SplitGeneSubset`, scatter `RunHiRes` over its File outputs, and add `MergeHiRes`. Merge by gene identity in original target-list order, using a disk-backed SQLite store to avoid loading all expression values into memory. Preserve numeric text and missing-value markers.

**Tech stack:** WDL 1.0, typed File inputs, Python standard library, the existing pinned official HiRes image, MiniWDL and Cromwell syntax validation, GitHub Actions.

**Spec:** The user's request in this chat and the design stated in commentary: `gene_subset` remains the full target list, `genes_per_chunk` defaults to 1000, each native call receives the full mixture, signature and fractions.

## Global constraints

- Target Terra managed Cromwell; WDL 1.0.
- No workflow-scope file-writing functions.
- Preserve incoming and generated File types until command rendering.
- Use named arguments and task-local newline-delimited File lists; generated list Files must stay File-typed.
- Check readable local Files and reject unresolved cloud URIs before computation.
- Keep logging, sample order, cell labels, and exact numeric text; preserve native 1 and NA markers.
- No local image build, authenticated analysis, or Terra submission.
- Default `genes_per_chunk=1000`; require genes_per_chunk > 0. Produce contiguous nonempty lists of at most that size; the last list can be smaller. The chunk count is ceil(target-gene count / genes_per_chunk).
- All native chunks receive identical full mixture, signature, fractions, and sampling/thread settings.
- Existing workflow output types stay unchanged; add chunk plans, lists, and per-chunk audit/log arrays.
- A completed synthetic smoke test does not establish native chunking equivalence or a successful Terra workflow run.

## Review focus

- Missing, repeated, unexpected, or overlapping genes must stop aggregation.
- Different samples, sample order, cell types, or window sizes must stop aggregation.
- Identical matrix basenames in different chunk directories must remain distinct localized Files.
- Gene and value text, including 1, NA, NaN and blank missing values, must be preserved.
- Newly generated task-output Files and command-time file lists must receive simulated cloud-to-local mapping in tests.

### Task 1: Split and merge interfaces

Files: `workflows/cell_type_specific_expression/cibersortx_hires.wdl`, `tests/cibersortx_hires/test_chunking.py`.

- [x] Add failing tests for fixed-size chunks, invalid sizes/lists, ordering, and generated File localization.
- [x] Add `SplitGeneSubset(File gene_subset, Int genes_per_chunk=1000)` with File array outputs and a structured plan.
- [x] Add `MergeHiRes` with original gene File, chunk-gene Files, matrix Files, input-validation Files and native-log Files. Reject inconsistent inputs and coverage; restore original target order per cell type.
- [x] Preserve all expression value text and use bounded memory during aggregation.
- [x] Run the new tests and existing native-task tests.

### Task 2: Workflow, documentation and cloud checks

Files: the HiRes WDL, existing HiRes tests and CI, examples, guide, workflow catalog, Cromwell test.

- [x] Wire split -> scatter -> merge. Forward unchanged full input Files into every native call.
- [x] Preserve existing File/Array[File] output types and expose per-chunk audit outputs.
- [x] Set `genes_per_chunk=1000` in the existing example; document all-gene lists and per-task resources.
- [x] Extend GitHub pinned-image tests to include split/merge and an end-to-end synthetic multi-chunk case.
- [x] Include the HiRes descriptor in Cromwell Womtool checks.
- [x] Run MiniWDL, file-scope/logging checks, task tests, image-routing/CI-selection tests, and review the final diff.
- [x] Push a reviewable branch/PR to prepare_QTL and verify the relevant GitHub checks; do not submit Terra.

## Verified result

- PR: https://github.com/AoU-Multiomics-Analysis/prepare_QTL/pull/87.
- All 26 local HiRes tests passed, including actual workflow-call expressions and simulated incoming/generated File localization.
- GitHub HiRes run 36966950298 passed host and pinned-image synthetic tests.
- GitHub descriptor run 36966950351 passed MiniWDL and Cromwell Womtool 87 checks.
- Image-release report passed; no image build was selected.
- Independent review found no remaining actionable issue.
- The chunked workflow and native equivalence remain untested on Terra.
