# Fast image releases

Use one PR for pipeline scripts, WDLs, documentation, and their existing tests.
The `release-ready` label starts image release when automation is enabled.
Only affected images are built or reused, and only affected WDL digests change.
Test-only changes do not rebuild images. Runtime tests do not run during release.

WDL validation stays automatic, including syntax, imports, logging, and Terra
file-scope checks. Image source fingerprints, registry digests, stage routing,
allowed pin paths, and stale-commit checks still protect the release.

To run optional tests, select a branch under GitHub Actions and use **Run workflow**:

- **Source unit checks**: source code checks using existing dependency images.
- **Cell-Type-Specific Expression CI**: cell-type analysis checks.
- **Pinned Image Smoke**: compact or full integration using that branch's pinned images.
- **Image Release Plan (Report Only)**: release-controller test suites on manual runs.

For tests of a new published image, select the PR branch after its digest update.
Tests can change in the same PR as scripts. No prerequisite test PR is needed.
Changes to privileged release policy under `ci/` or `.github/` still merge separately.
Tests do not block publishing or pin commits. A published image can be untested;
run manual checks when the change warrants them.

The dispatcher ignores verified digest-only App commits and skips superseded
release requests. There is no automatic PR merge. No complete Terra run is implied
by publishing an image or passing local or GitHub checks.
