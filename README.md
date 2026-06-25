# prepare_QTL

Pipeline for preparing molecular phenotype and genotype data for QTL (Quantitative Trait Loci) analysis. Supports expression (eQTL), splicing (sQTL), and proteomics (pQTL) data types.

---

## Repository Structure

```
prepare_QTL/
├── scripts/        # R scripts for data processing and normalization
├── workflows/      # WDL workflows for running analyses on a cloud platform
└── envs/           # Docker environments
```

---

## Scripts

All scripts are written in R and are invoked via the command line using `Rscript`. They are also bundled in the Docker image defined in `envs/PhenotypePCs/Dockerfile`.

### `scripts/PrepareExpression.R`

Prepares RNA-seq gene expression data for eQTL analysis.

**What it does:**
1. Loads raw count data from a GCT-formatted file (`.gct` or tab-separated).
2. Loads a GENCODE GTF annotation file and extracts transcription start site (TSS) positions for each gene.
3. Filters to a specified list of samples.
4. Filters genes by expression count (retains genes with counts > 6 in at least 20% of samples).
5. Normalizes counts using edgeR TMM normalization followed by CPM transformation.
6. Applies rank-based inverse normal (RankNorm) transformation to each gene by default, or centers and scales each gene when rank normalization is disabled.
7. Merges the normalized expression values with TSS locations and writes the result to a compressed BED file (`.expression.bed.gz`).

**Inputs:**
- `--CountGCT`: GCT or TSV file of raw RNA-seq count data.
- `--AnnotationGTF`: GENCODE GTF file used to extract TSS locations.
- `--SampleList`: File containing the list of sample IDs to include.
- `--OutputPrefix`: Prefix for the output file.
- `--RankNormalize`: Whether to apply RankNorm transformation (`true`, default) or only center and scale (`false`).

**Output:** `<OutputPrefix>.expression.bed.gz`

---

### `scripts/PrepareProteomics.R`

Prepares Olink proteomics data for pQTL analysis.

**What it does:**
1. Reads Olink proteomics data from a TSV, gzipped TSV, or Parquet file.
2. Loads a GENCODE GTF file to extract TSS positions.
3. Uses Ensembl BioMart to map UniProt IDs to Ensembl gene IDs.
4. Filters to a specified list of samples.
5. Pivots data to wide format (proteins as columns, samples as rows).
6. Applies rank-based inverse normal (RankNorm) transformation to each protein by default, or centers and scales each protein when rank normalization is disabled.
7. Merges with TSS locations and writes the result to a compressed BED file (`.protein.bed.gz`).

**Inputs:**
- `--ProteomicData`: TSV, TSV.gz, or Parquet file of normalized Olink protein expression data.
- `--AnnotationGTF`: GENCODE GTF file used to extract TSS locations.
- `--SampleList`: File containing the list of sample IDs to include.
- `--OutputPrefix`: Prefix for the output file.
- `--RankNormalize`: Whether to apply RankNorm transformation (`true`, default) or only center and scale (`false`).

**Output:** `<OutputPrefix>.protein.bed.gz`

---

### `scripts/PrepareSpliceData.R`

Prepares LeafCutter splice junction data for sQTL analysis.

**What it does:**
1. Loads a BED file of splice junction quantifications produced by LeafCutter.
2. Filters columns to a specified list of samples.
3. Applies rank-based inverse normal (RankNorm) transformation to each splice junction by default, or centers and scales each splice junction when rank normalization is disabled.
4. Sorts junctions by chromosome, start, end, and Ensembl gene ID.
5. Writes the normalized splice data to a compressed BED file (`.splicing.bed.gz`).

**Inputs:**
- `--SpliceData`: BED file of splice junction quantifications from LeafCutter.
- `--SampleList`: File containing the list of sample IDs to include.
- `--OutputPrefix`: Prefix for the output file.
- `--RankNormalize`: Whether to apply RankNorm transformation (`true`, default) or only center and scale (`false`).

**Output:** `<OutputPrefix>.splicing.bed.gz`

---

### `scripts/calculate_PCs.R`

Computes phenotype principal components (PCs) from a normalized BED file.

**What it does:**
1. Reads a normalized BED file (e.g., expression, splicing, or proteomics).
2. Runs PCA on the molecular phenotype data using `PCAtools::pca`.
3. Selects the optimal number of PCs using the Gavish-Donoho method (`chooseGavishDonoho`).
4. Writes the phenotype PCs to a TSV file.

**Inputs:**
- `--bed_file`: A normalized BED file of molecular phenotype data.
- `--output_prefix`: Prefix for the output file.

**Output:** `<output_prefix>_phenotype_PCs.tsv`

---

### `scripts/MergeCovariates.R`

Merges genotype PCs and molecular phenotype PCs into a single covariate file for QTL analysis.

**What it does:**
1. Reads genotype PCs and molecular phenotype PCs from TSV files.
2. Inner-joins the two tables by sample ID.
3. Transposes the merged table so that rows are covariate names and columns are sample IDs, as required by tensorQTL.
4. Writes the combined covariate matrix to a TSV file.

**Inputs:**
- `--GenotypePCs`: TSV file of genotype principal components.
- `--MolecularPCs`: TSV file of molecular (phenotype) principal components.
- `--OutputPrefix`: Prefix for the output file.

**Output:** `<OutputPrefix>_QTL_covariates.tsv`

---

## Workflows

All workflows are written in WDL (Workflow Description Language) and are designed to run on a cloud platform such as Terra. Each workflow wraps one or more tasks that call the scripts above or external tools.

### `workflows/prepare_eQTL.wdl`

End-to-end workflow for preparing gene expression data for eQTL analysis.

**Steps:**
1. Runs `PrepareExpression.R` to normalize expression data and produce a BED file.
2. Runs `calculate_PCs.R` (via `calculate_phenotypePCs.wdl`) to compute phenotype PCs from the expression BED file.
3. Optionally runs `MergeCovariates.wdl` when `AdditionalCovariates` is provided.

**Inputs:** Raw count GCT file, GENCODE GTF, sample list, output prefix, rank normalization toggle, optional additional covariates TSV, resource parameters.
**Outputs:** Expression BED file (`.expression.bed.gz`), phenotype PCs TSV, and optionally merged QTL covariates TSV.

---

### `workflows/prepare_pQTL.wdl`

End-to-end workflow for preparing Olink proteomics data for pQTL analysis.

**Steps:**
1. Runs `PrepareProteomics.R` to normalize protein expression data and produce a BED file.
2. Runs `calculate_PCs.R` (via `calculate_phenotypePCs.wdl`) to compute phenotype PCs from the protein BED file.
3. Optionally runs `MergeCovariates.wdl` when `AdditionalCovariates` is provided.

**Inputs:** Olink proteomics data file, GENCODE GTF, sample list, output prefix, rank normalization toggle, optional additional covariates TSV, resource parameters.
**Outputs:** Protein BED file (`.protein.bed.gz`), phenotype PCs TSV, and optionally merged QTL covariates TSV.

---

### `workflows/prepare_sQTL.wdl`

End-to-end workflow for preparing splice junction data for sQTL analysis.

**Steps:**
1. Runs `PrepareSpliceData.R` to normalize LeafCutter splice data and produce a BED file.
2. Runs `calculate_PCs.R` (via `calculate_phenotypePCs.wdl`) to compute phenotype PCs from the splicing BED file.
3. Optionally runs `MergeCovariates.wdl` when `AdditionalCovariates` is provided.

**Inputs:** LeafCutter BED file, sample list, output prefix, rank normalization toggle, optional additional covariates TSV, resource parameters.
**Outputs:** Splicing BED file (`.splicing.bed.gz`), phenotype PCs TSV, and optionally merged QTL covariates TSV.

---

### `workflows/calculate_phenotypePCs.wdl`

Workflow that computes phenotype PCs from any normalized molecular phenotype BED file.

**Steps:**
1. Runs `calculate_PCs.R` on the input BED file using the Gavish-Donoho method to select the number of PCs.

**Inputs:** Normalized BED file, output prefix, resource parameters.  
**Outputs:** Phenotype PCs TSV (`<OutputPrefix>_phenotype_PCs.tsv`).

---

### `workflows/MergeCovariates.wdl`

Workflow that merges additional covariates, such as genotype PCs, and molecular phenotype PCs into a single covariate file ready for tensorQTL.

**Steps:**
1. Runs `MergeCovariates.R` to inner-join genotype and molecular PCs and transpose the result.

**Inputs:** Additional covariates TSV with `sample_id`, molecular PCs TSV, output prefix.
**Outputs:** Combined QTL covariate file (`<OutputPrefix>_QTL_covariates.tsv`).

---

### `workflows/prepare_VCF.wdl`

Comprehensive genotype preparation workflow for All of Us (AoU) data starting from a Hail MatrixTable.

**Steps:**
1. Imports the `FilterMTAndExportToVCF` workflow from the [MTtoVCF](https://github.com/AoU-Multiomics-Analysis/MTtoVCF) repository to filter a Hail MatrixTable by sample list and allele count threshold and export a VCF.
2. Converts the VCF to PLINK2 binary format (`.pgen`/`.pvar`/`.psam`) using PLINK2.
3. Computes genotype PCs from the VCF using an R script.

**Inputs:** Hail MatrixTable URI, sample list, allele count threshold, maximum allele length, output prefix, cloud paths, genotype PCs R script.  
**Outputs:** VCF file, PLINK2 files (`.pgen`, `.pvar`, `.psam`), genotype PCs TSV.

---

### `workflows/convertVCF2Plink.wdl`

Workflow that converts a VCF file to PLINK2 binary format.

**Steps:**
1. Runs `plink2 --make-pgen` to convert a VCF to PLINK2 `.pgen`/`.pvar`/`.psam` format, setting variant IDs to `CHR:POS_REF_ALT` and filtering to autosomes (chr1–22).

**Inputs:** VCF file, output prefix, maximum allele ID length.  
**Outputs:** PLINK2 files (`.pgen`, `.pvar`, `.psam`).

---

### `workflows/calculateGenotypePCs.wdl`

Workflow that computes genotype principal components from a VCF file.

**Steps:**
1. Runs an R script (supplied as input) against the VCF to compute genotype PCs.

**Inputs:** VCF file, output prefix, genotype PCs R script.  
**Outputs:** Genotype PCs TSV (`<OutputPrefix>_genetic_PCs.tsv`).

---

### `workflows/calculateAF.wdl`

Workflow that calculates allele frequencies from PLINK2 genotype files for a specified set of samples.

**Steps:**
1. Optionally strips a header row from the sample list file.
2. Runs `plink2 --freq` on the provided PLINK2 files restricted to the given sample list.

**Inputs:** PLINK2 files (`.pgen`, `.pvar`, `.psam`), sample list, output prefix, resource parameters.  
**Outputs:** Allele frequency file (`<prefix>.afreq`).

---

### `workflows/calculateGenotypeDosage.wdl`

Workflow that extracts genotype dosage values from a VCF file.

**Steps:**
1. Uses `bcftools query` to extract sample IDs and `bcftools +dosage` to compute dosage from genotype (GT) fields.
2. Prepends a header row (CHROM, POS, REF, ALT, sample IDs) and compresses the output with bgzip.
3. Indexes the output with tabix.

**Inputs:** VCF file (`.vcf.gz`), number of threads.  
**Outputs:** Bgzipped dosage matrix (`<vcf_basename>.dose.tsv.gz`) and its tabix index.

---

## Docker Environment

The `envs/PhenotypePCs/Dockerfile` defines the Docker image used by most WDL tasks (published as `ghcr.io/aou-multiomics-analysis/prepare_qtl:main`). It is built automatically on every push or pull request to `main` via the GitHub Actions workflow in `.github/workflows/docker-image.yml`.

The image includes the following R packages:
- `tidyverse`, `data.table`, `arrow`, `optparse`, `janitor`
- `PCAtools`, `RNOmni`, `edgeR`
- `biomaRt`, `biomaRtr`, `rtracklayer`
- `patchwork`
