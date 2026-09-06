# Per-cell-type BED filtering

`CellTypeDeconvolution` scatters `FilterCellTypeBed` over the actual exported
`Array[File]`. Each job localizes one BED, the small export inventory, and the
optional prepared Haemopedia summary. Haemopedia is still prepared once.

The job selects exactly one inventory row by BED basename, then checks the BED
checksum and dimensions. It applies the existing negative-value, expression,
and optional residual rules. Expression remains in linear CPM. A gene with a
negative value is removed only from the affected cell type.

`MergeFilterReports` runs after all shards finish. It restores export-inventory
order, checks complete cell-type coverage and identical sample order, combines
the tables, and recreates the cross-cell-type negative plots and reference plots.
It reads only reports and sample lists, not the full BED matrices. BED files stay
at their shard output locations. Public output names and prepare-eQTL wiring are
unchanged. The existing `FilterCellTypeBeds` call label now identifies this merge
step; `FilterCellTypeBed` identifies the per-cell-type work.

Files remain typed as `File` or `Array[File]`. Required report paths are written
into task-local lists at command rendering, after localization. Inventory paths
are metadata basenames and are never used to reconstruct cloud input paths.

Existing memory and disk settings apply to each filtering shard. The merge uses
the existing summary memory setting and a 20 GB disk. Parallel jobs can reduce
elapsed time but increase simultaneous VM and disk use. A failed shard fails the
workflow; no cell type is silently omitted. Normal Cromwell cache rules still
apply, including the shared inventory input; scattering does not fix cache misses.

Local validation compares serial and scattered results with no reference, with
reference filtering, and with residual filtering. It also checks reversed shard
order, missing shards, sample-order mismatches, CLI execution, and simulated
cloud-to-local file mapping. No complete Terra run has been performed for this
change. No local Docker build is required.
