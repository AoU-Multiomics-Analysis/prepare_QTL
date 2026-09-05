# Directory-based script stages

## Goal

Adding a script inside a registered stage directory must select that stage for
image release without another per-file registry entry. A new WDL task must use
the image input for the stage of the script it runs. Preserve scientific
calculations, file formats, WDL 1.0, and Terra file localization.

## Approach

Use directory rules, explicit WDL image inputs, and static validation. Do not
generate WDL tasks or infer dependencies from arbitrary shell or R code.

Keeping per-file registration would preserve the present maintenance burden.
Generating all WDL wiring from a new framework would add unnecessary scope.
Directory rules fit the existing release planner and image registry.

## Layout and scope

Keep existing stage IDs, image repositories, and image-input names. Organize
cell-type scripts under these directories:

```text
scripts/cell_type_specific_expression/
  estimation/
  fit/
  export/
  downstream/
  shared/
  tools/
```

Put a stage's entry scripts and private R modules in its directory. Put modules
used by more than one stage in shared; initially changes there select all four
cell-type stages. This is conservative and avoids missing a dependent stage.
Preserve the current module load order rather than relying on a new recursive
alphabetical order. Keep development-only tools separate from runtime code.

Other existing single-stage directories, such as proteomics and splicing,
remain in place. Separate RNA-SeQC scripts from expression preparation under
the expression family so both stages can use directory rules. Keep Rust crate
directories as their stage roots. Shared common scripts keep their explicit
list of consuming stages. This is a repository-wide registry convention, not
a requirement to rename directories that already form a stage boundary.

## Stage images and shared bases

The user approved separate stage images on 2026-09-05. Each stage will have its
own image repository and immutable digest default. Stages with compatible
dependencies share a published base image pinned by digest. Bases contain
dependencies, not the repository's analysis scripts. A stage layer contains
its own code, its required shared modules, and compatibility entry points.
Different languages or incompatible dependency sets use different bases.

Do not rebuild the dependency base for a stage-script change. A dependency or
base recipe change selects every child stage. Include the resolved base digest
in each child build fingerprint; the same source with a different base is a
different image. Never use a mutable base tag for a released child image.
Resolve and verify base images before child builds. Retain both base and child
digests. The actual publication needs release-publish approval.

Stage-image repositories, child/base build ordering, and pin-repository
transition rules must be installed as trusted policy before the source/image
migration release. The first implementation plan covers the independently
testable directory and routing policy. A second plan must cover this build
graph and its security review before the image migration begins. It must not
be presented as implemented by the first plan.

## Registration and validation

Use directory globs in ci/image-stages.yml. New paths under a registered root
inherit its stage. Build-path coverage remains separate from stage selection.
Keep ci/release-pins.yml as the explicit list of image-default edit targets.

Validate that runtime scripts belong to a registered stage or shared root,
that roots do not assign conflicting ownership, and that selected stage files
are covered by their image's build inputs. Shared consumers and development
exclusions remain explicit. New stages still require registry and pin entries.

Extend existing MiniWDL-based image-routing tests. For supported literal
repository script invocations, resolve the script's stage from its path and
check the image passed through workflow calls to the task runtime. Do not use
a growing table of task names for these tasks. Detect incorrect forwarding,
hard-coded task images that bypass the stage input, and unknown script roots.

Keep explicit validation for existing inline-only tasks and external images.
An unsupported dynamic script path must fail with instructions to add a
reviewed contract; it must not be silently treated as correctly wired. Static
validation does not prove arbitrary imported module dependencies or scientific
correctness. New tasks still require runtime tests.

## Compatibility and rollout

Do not change existing WDL command paths merely to match source directories.
New images must provide the old entry-point paths through thin compatibility
entry points or build-time layout mapping, with one canonical implementation.
The implementation plan must choose the mapping per image and test both old
entry points and new stage paths. Do not duplicate scientific implementations.

This keeps current main WDLs usable with their current pinned images. New task
commands that require new paths cannot be used with an old digest; add them
only with a tested candidate image. Do not change existing digests without
publication and runtime tests. Path and Dockerfile changes alter fingerprints
and can require rebuilding multiple images once during migration.

Use two release phases because trusted release policy comes from main:

1. Add directory rules, validators, compatibility contracts, and tests while
   still accepting the current registered paths. Merge this policy separately.
2. Move source files and update container mappings in a source migration PR.
   Use the merged policy to build and test candidate images and update pins.

Remove temporary legacy source classifications after migration in a small
policy cleanup. Runtime compatibility paths can remain; they are not a second
source implementation. Never relax the release controller's restriction on
candidate PRs changing their own trusted tests or policy.

## Tests and acceptance

- A previously unknown script under each stage root selects the correct stage.
- Shared changes select all declared consumers; outside-root scripts fail.
- A new synthetic WDL task using its correct stage image passes without adding
  its task name to a routing table. A wrong image or unresolved command fails.
- Existing parent/child image inputs and immutable defaults remain valid.
- Compatibility entry points load the same implementation and preserve CLI
  arguments, module order, outputs, and logging.
- Existing Python, R, WDL, and image-release checks still run. Add container
  checks to GitHub Actions; do not build Docker images locally.

Implementation does not enable releases, publish images, push branches, merge
PRs, or submit Terra jobs without the applicable user approval. Complete Terra
execution and call-cache behavior remain separate, untested claims.
