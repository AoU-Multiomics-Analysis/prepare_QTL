# scDblFinder across the complete reference cohort

Purpose: score all eligible capture libraries, exclude predicted doublets, and recreate the CIBERSORTx single-cell input reference.

The scope preserves the existing reference rules: 22 capture libraries from 8 donors; 10x 3-prime v3 only; TSP10 remains excluded. Detection uses all 71,101 available cells in those captures, including 643 context cells outside the reviewed 70,458-cell reference. The context cells are not added to the output reference. Each capture is scored separately using original raw counts from the Tabula Sapiens Blood H5AD. All-zero features are dropped for scoring; remaining raw values are unchanged.

Use the same Bioconductor scDblFinder 1.14.0 installation and R 4.3.1 as the TSP7 pilot. The R detection script is unchanged. Use clusters=FALSE, the automatic expected doublet rate, package defaults, and a fixed seed for each capture. Preserve the pilot's TSP7 seeds. Store scores, calls, method metadata, and package versions. No additional low-quality-cell filter is applied. SoupX and CellBender are not run.

Recreate the reference by removing only predicted doublet columns from the existing reference with sex chromosome genes excluded. Preserve reviewed labels, cell order, integer count values, bulk-detected gene selection, the prior 17 gene exclusions, and the 741 sex chromosome gene exclusions. The 252 NKT review cells and one unresolved review cell remain excluded. Drop genes that become zero across all retained reference cells, and report them. The original reference remains unchanged. Export gzip directly, with a selected-cell manifest, removed-cell manifest, validation, and file hashes.

The config contains absolute existing local input paths, the capture plan, seed mapping, and named Python/R executable and package-library paths. Scripts receive named CLI arguments. No Docker image or cloud executor is used. Snakemake runs with four cores; each detection task requests four cores, so captures are processed one at a time. The engine stages this source; inputs, package libraries, counts, and results remain outside it.

Input readiness requires the original H5AD, current cell manifest, flagged 226 CD8 cell manifest, focused marker table, existing sex-filtered reference and its validation record, host Python with numpy/scipy/h5py/pandas, and the existing private R library with scDblFinder 1.14.0. Test scripts check that reference call alignment rejects missing or duplicate cells and that reference filtering preserves labels, order, and count values.

The endpoint is a new reference input, not a newly derived signature matrix. Marker derivation and S-mode fraction estimation are not run. This workflow is local; it has not been run on Terra. Predicted doublets are model calls, not independent proof of poor cell quality, and residual contamination or annotation errors can remain.
