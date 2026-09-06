# Fast image releases

Script changes still trigger image deployment for a same-repository PR with the
`release-ready` label, when release automation is enabled. The release builds
only the image families required by affected stages. It tests the images and
updates only the affected WDL digest defaults. Unaffected defaults stay pinned.

Automatic cell-type releases run one small test file per affected stage:

| Stage | Test |
| --- | --- |
| Estimation | HSPE batches |
| Fit | TCA fit, cleanup, and export alignment |
| Export | BED outputs |
| Downstream | Reference filtering |

Other image families keep their existing short task tests. Dependency changes
select their consumer stages, so those tests also check that the image can start
and load its packages. They do not start a full cell-type workflow.

WDL changes keep syntax, import, routing and Terra file-scope checks. A WDL-only
change does not build an image or select runtime tests for unrelated stages.

Broad source tests and compact/full integration tests are manual GitHub Actions:
`Source unit checks`, `Cell-Type-Specific Expression CI`, and `Pinned Image Smoke`.
Use integration tests after changes to connections between tasks, or when needed
to investigate a failure. These tests remain available; they are not release gates.

The dispatcher ignores the release App's verified digest-only commit. Each
automatic request also records its trigger commit. A queued request stops before
building or testing if the PR head has changed. Manual release requests can omit
the expected commit to retry the current head.

These checks do not prove that a complete workflow runs on Terra. Pinned image
digests preserve the selected software; retain workflow revisions and run inputs
as well to reproduce an analysis.
