# Release pin proposal implementation plan

**Goal:** Produce a reviewable patch for affected stage defaults without publishing images, changing working WDLs, or writing to GitHub.

**Architecture:** A separate target registry maps release stages to WDL input defaults. The proposal tool reads committed sources at exact base/head revisions, calls the dependency planner, validates candidate image references, and edits only AST-located String literals in memory. It writes a patch and a release record into a new output directory.

**Spec:** `docs/stage-image-pinning-draft.md`, first part of the build/test/pin sequence. This artifact is a candidate, not permission to merge or evidence of runtime validation.

## Tasks

- [ ] Add `tests/test_release_pin_proposals.py` with synthetic WDLs. Check downstream-only updates, unchanged fitting defaults, duplicate/missing targets, invalid repositories/tags, inconsistent defaults, unsupported stages, and no-op behavior.
- [ ] Run `python tests/test_release_pin_proposals.py` and confirm failures before implementation.
- [ ] Add `ci/release-pins.yml` for the five stage inputs currently exposed by the two cell-type entry points. Other stage families must fail explicitly if a proposal tries to update them.
- [ ] Add `ci/propose_image_pins.py`: `propose(config, targets, plan, files, candidates)` returns changed file text only. Use MiniWDL AST source positions, not a global text replacement. Validate all target defaults before proposing any edit.
- [ ] Add CLI checks against exact expected base/head SHA values. Read configuration and WDLs from the selected committed head, not uncommitted files. Require base to be an ancestor of head. Write only to a newly created output directory.
- [ ] Test the CLI against a temporary Git repository: accepted proposal yields a patch, unchanged files stay unchanged, stale refs yield no artifacts, and the tool never changes tracked files.
- [ ] Add tests to read-only CI and document the CLI, trust boundary, limitations, and remaining build/test/updater integration.
- [ ] Run tests, actionlint, and diff checks. Request code review. Do not enable a GitHub write job, create credentials, or submit Terra jobs.

## Release boundary

The tool validates the syntax and repository of supplied digests. It does not certify registry existence, image source provenance, smoke-test success, or GitHub event authenticity. Those checks belong in the future build/test and trusted-updater jobs. Local ref checks must not be described as proof that a remote PR is current.

## Implementation checkpoint

Implemented the candidate patch tool, target registry, CLI, read-only CI tests,
and user documentation on `codex/release-pin-proposals`.

Seventeen proposal tests pass, including temporary Git repository tests for
stale head/base refs, unchanged tracked files, committed configuration, and
patch application. Eighteen existing image/planner/CI-selection tests also pass.
Actionlint and diff checks pass. Review found the missing-final-newline patch
edge case; it was fixed and covered by a regression test.

No build/test orchestration or trusted commit job is enabled by this checkpoint.
No production WDL was changed by this layer, and no GitHub or Terra write was
performed. The build/test/commit integration remains the next implementation step.
