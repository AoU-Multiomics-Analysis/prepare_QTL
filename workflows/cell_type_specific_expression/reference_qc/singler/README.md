# SingleR review of retained Tabula Sapiens blood labels

This local Snakemake workflow reviews all 65,191 cells retained after scDblFinder. It does not change labels or rebuild the CIBERSORTx input.

SingleR 2.2.0 and celldex 1.10.1 run under R 4.3.1. Monaco Immune main and fine labels and Human Primary Cell Atlas main labels are classified independently. Raw UMI counts are valid for SingleR rank correlations. Each capture is read as a sparse matrix. Duplicate gene symbols are summed. Features associated with X or Y in GENCODE v48 are excluded. Globin genes are kept for erythroid label review. Training uses default classic marker selection, shared gene symbols, exact neighbors, four threads, and default classification settings. Score pruning is applied once per model across all retained cells.

The typed file paths in the configuration refer to local, verified raw-count capture exports from the completed scDblFinder run. The feature map and two downloaded reference RDS files are external inputs. Configuration is local Snakemake configuration, not a WDL argument wrapper. This workflow has not run on Terra.

Inputs: counts_dir with counts.mtx.gz/features.tsv/cells.tsv per capture; capture_plan; selected retained cell table; feature_map with feature_id/feature_symbol/exclude_sex; reference_monaco; reference_hpca; flagged original CD8 review table; CLC_review six-cell marker table; rscript, r_library, python, captures. Scripts use named CLI arguments and reject missing inputs or cell mismatches.

Outputs: per-cell raw/pruned calls, delta scores, full score matrices and prediction RDS files; raw marker panel; capture/training metadata; disagreement and uncertainty lists; focused six CLC-high and 51 retained S100-high cell reviews; status for all original 226 flagged CD8 cells; by-donor count tables; session information.

Coverage limits: Monaco has no erythroid or platelet reference. Its main B label covers plasmablasts; only its fine labels assess plasma compatibility. HPCA T labels assess the broad T lineage, not CD4/CD8 identity. MAIT/gamma-delta calls do not resolve the CD4/CD8 subtype. A pruned call is uncertain. None of these calls is a probability or an automatic instruction to remove a cell.

Tests: python -m unittest discover -s tests. Execute scientific work only through the approved NGS Workbench plan.
