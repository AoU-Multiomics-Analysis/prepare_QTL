# Image release planner: report-only trial

This is the first part of the stage-pinning design. It reports affected release
groups across the repository. It does not pin images, change WDLs, build containers,
publish packages, open PRs, or submit Terra jobs.

## Try a filtering change

```sh
python ci/plan_image_updates.py \
  --changed scripts/cell_type_specific_expression/filter_cell_type_beds.R
```

Expected: image build dependencies are `cell_type` and `standard`; the only
proposed pin-update group is `cell_downstream`.

Why two builds? The cell-type image copies its script directory and the standard
image copies all scripts. Rebuilding an image does not mean every stage must
adopt that new image. The existing build workflows remain unchanged in this trial.

For a real branch comparison:

```sh
python ci/plan_image_updates.py --base origin/main --head HEAD
python ci/plan_image_updates.py --validate
```

The new read-only GitHub workflow publishes the report in its job summary for
each PR. No bot token or secrets are needed for this report. A failed classification
check does not change branch protection unless an administrator makes it required.

## Repo-wide groups

The map includes cell-type estimation, fitting, export, and downstream processing;
standard expression, proteomics, splicing, and methylation; Rust methylation;
and RNA-SeQC aggregation. Common QTL helpers affect all standard QTL groups.

Genotype workflows use third-party images. Their WDL changes are reported for
review; this system does not build or automatically update those external images.
Their digest selection remains a separate, explicit operation.

These are planned release groups, not newly installed WDL inputs. Some existing
workflows still use shared or hard-coded image references. They need a later
interface update before the groups can control their runtime images independently.

## Update the map

`ci/image-stages.yml` defines image build paths, stage source paths, shared helpers,
explicitly ignored development tools, and workflow families. Patterns use Python
fnmatch semantics: an asterisk may match path separators.

Unknown source paths under scripts, rust, or envs fail classification. New WDL
families also fail. Add their actual consumers and a regression case before
accepting them. A new file within an existing broad family rule inherits that
family; this does not prove its fine-grained dependencies have been audited.

WDL-only changes request validation, not image builds or pin updates. Test-only
changes are listed without proposing image changes. Renames are treated as a
deletion and addition so both paths can affect selection.

The current map is a tested proposal, not automatic R dependency
analysis. Before enabling pin updates, audit shared-function callers and test the
mixed image set. The bootstrap currently loads every R module, which requires
particular care with top-level side effects and changes to shared functions.

## What remains before automatic pinning

1. Audit and approve the dependency map using the reports from actual PRs.
2. Select verified initial image digests and add the stage inputs to WDLs.
3. Test the actual mixture of unchanged and new images in GitHub Actions.
4. Add a trusted pin updater with stale-result checks and no-op detection.
5. Approve the GitHub App installation, required checks, and image retention.
6. With user approval, verify caching across otherwise identical Terra runs.

Keep using the WDL source on main. Once the later updater exists, it will commit
stage digest changes in the code PR. It will not require a new WDL URL. Existing
Terra submissions and explicit input overrides still keep their own settings.
