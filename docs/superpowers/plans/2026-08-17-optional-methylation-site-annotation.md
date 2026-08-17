# Optional methylation site annotation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make methylation site annotation optional while keeping it enabled by default and increase its default resources to 256 GB memory and 200 GB disk.

**Architecture:** Gate the existing annotation call with an `AnnotateSites` Boolean in the public cohort workflows. Propagate the Boolean and resource defaults through the manifest-based wrapper and array-based workflow; expose an optional annotation output because disabled runs do not produce an annotation file.

**Tech Stack:** WDL 1.0, miniwdl, shell regression tests.

## Global Constraints

- `AnnotateSites` defaults to `true`.
- `AnnotationMemoryGB` defaults to `256`.
- `AnnotationDiskGB` defaults to `200`.
- Annotation remains downstream-independent; disabling it must not remove aggregation, connectivity, PCs, or QTL-covariate outputs.
- Do not alter the annotation task command or output filename.

---

### Task 1: Add failing workflow regression coverage

**Files:**
- Create or modify: `tests/methylation/test_annotation_options.sh`

**Interfaces:**
- Consumes: `workflows/methylation/AggregateMethylationCohort.wdl`, `workflows/methylation/AggregateMethylationCohortArrays.wdl`, and `workflows/methylation/merge_methylation.wdl`.
- Produces: assertions for the default-enabled Boolean, 256/200 resource defaults, conditional annotation call, and optional annotation output.

- [ ] **Step 1: Write the failing test**

Use `miniwdl check` and text assertions to require `AnnotateSites = true`, `AnnotationMemoryGB = 256`, `AnnotationDiskGB = 200`, an annotation `if` block, and `File? PassingSiteAnnotations` in `cohort_aggregation.wdl`, `AggregateMethylationCohortArrays.wdl`, `AggregateMethylationCohort.wdl`, and `merge_methylation.wdl`.

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
bash tests/methylation/test_annotation_options.sh
```

Expected: failure because the current workflows have no `AnnotateSites` input, use 64/100 defaults, and call annotation unconditionally.

- [ ] **Step 3: Commit the failing test**

```bash
git add tests/methylation/test_annotation_options.sh
git commit -m "test: define optional methylation annotation behavior"
```

### Task 2: Implement optional annotation and resource defaults

**Files:**
- Modify: `workflows/methylation/AggregateMethylationCohortArrays.wdl`
- Modify: `workflows/methylation/AggregateMethylationCohort.wdl`
- Modify: `workflows/methylation/cohort_aggregation.wdl`
- Modify: `workflows/methylation/merge_methylation.wdl`
- Test: `tests/methylation/test_annotation_options.sh`

**Interfaces:**
- Public inputs: `Boolean AnnotateSites = true`, `Int AnnotationMemoryGB = 256`, and `Int AnnotationDiskGB = 200`.
- Internal workflow input: pass all three values through to the existing annotation call.
- Public output: `File? PassingSiteAnnotations` for workflows where the call is conditional.

- [ ] **Step 1: Add the conditional call in the array workflow**

Add `AnnotateSites` to the array workflow, wrap its existing `AnnotateMethylationSites` call in `if (AnnotateSites) { ... }`, pass the resource inputs unchanged, and expose `File? PassingSiteAnnotations = AnnotateMethylationSites.PassingSiteAnnotations`.

- [ ] **Step 2: Add the conditional call in the manifest aggregation workflow**

Add `Boolean AnnotateSites = true` to `AggregateMethylationData`, wrap `Annotation.AnnotateMethylationCohortSites as AnnotateSites` in `if (AnnotateSites)`, and expose `File? PassingSiteAnnotations`.

- [ ] **Step 3: Propagate the option through public wrappers**

Add and forward `AnnotateSites`, `AnnotationMemoryGB = 256`, and `AnnotationDiskGB = 200` through `AggregateMethylationCohort.wdl` and `merge_methylation.wdl`, and make their annotation outputs optional.

- [ ] **Step 4: Run the focused test**

Run:

```bash
bash tests/methylation/test_annotation_options.sh
```

Expected: PASS for all three workflows.

- [ ] **Step 5: Commit the implementation**

```bash
git add workflows/methylation/AggregateMethylationCohortArrays.wdl workflows/methylation/AggregateMethylationCohort.wdl workflows/methylation/cohort_aggregation.wdl workflows/methylation/merge_methylation.wdl tests/methylation/test_annotation_options.sh
git commit -m "feat: make methylation site annotation optional"
```

### Task 3: Validate WDL interfaces and repository cleanliness

**Files:**
- No additional source changes expected.

- [ ] **Step 1: Validate all affected descriptors**

```bash
miniwdl check workflows/methylation/AggregateMethylationCohortArrays.wdl
miniwdl check workflows/methylation/AggregateMethylationCohort.wdl
miniwdl check workflows/methylation/merge_methylation.wdl
git diff --check
```

- [ ] **Step 2: Confirm only intended changes are present**

```bash
git status --short
git diff --stat HEAD~2..HEAD
```
