# Genotype Workflows

[Back to main README](../README.md)

These WDL workflows prepare or summarize genotype data for downstream QTL analysis.

## `workflows/genotype/prepare_VCF.wdl`

Comprehensive genotype preparation workflow for All of Us (AoU) data starting from a Hail MatrixTable.

**Steps:**
1. Imports the `FilterMTAndExportToVCF` workflow from the [MTtoVCF](https://github.com/AoU-Multiomics-Analysis/MTtoVCF) repository to filter a Hail MatrixTable by sample list and allele count threshold and export a VCF.
2. Converts the VCF to PLINK2 binary format (`.pgen`/`.pvar`/`.psam`) using PLINK2.
3. Computes genotype PCs from the VCF using an R script.

**Inputs:** Hail MatrixTable URI, sample list, allele count threshold, maximum allele length, output prefix, cloud paths, genotype PCs R script.

**Outputs:** VCF file, PLINK2 files (`.pgen`, `.pvar`, `.psam`), genotype PCs TSV.

## `workflows/genotype/convertVCF2Plink.wdl`

Workflow that converts a VCF file to PLINK2 binary format.

**Steps:**
1. Runs `plink2 --make-pgen` to convert a VCF to PLINK2 `.pgen`/`.pvar`/`.psam` format, setting variant IDs to `CHR:POS_REF_ALT` and filtering to autosomes (`chr1` through `chr22`).

**Inputs:** VCF file, output prefix, maximum allele ID length.

**Outputs:** PLINK2 files (`.pgen`, `.pvar`, `.psam`).

## `workflows/genotype/calculateGenotypePCs.wdl`

Workflow that computes genotype principal components from a VCF file.

**Steps:**
1. Runs an R script, supplied as input, against the VCF to compute genotype PCs.

**Inputs:** VCF file, output prefix, genotype PCs R script.

**Outputs:** Genotype PCs TSV (`<OutputPrefix>_genetic_PCs.tsv`).

## `workflows/genotype/calculateAF.wdl`

Workflow that calculates allele frequencies from PLINK2 genotype files for a specified set of samples.

**Steps:**
1. Optionally strips a header row from the sample list file.
2. Runs `plink2 --freq` on the provided PLINK2 files restricted to the given sample list.

**Inputs:** PLINK2 files (`.pgen`, `.pvar`, `.psam`), sample list, output prefix, resource parameters.

**Outputs:** Allele frequency file (`<prefix>.afreq`).

## `workflows/genotype/calculateGenotypeDosage.wdl`

Workflow that extracts genotype dosage values from a VCF file.

**Steps:**
1. Uses `bcftools query` to extract sample IDs and `bcftools +dosage` to compute dosage from genotype (`GT`) fields.
2. Prepends a header row (`CHROM`, `POS`, `REF`, `ALT`, sample IDs) and compresses the output with bgzip.
3. Indexes the output with tabix.

**Inputs:** VCF file (`.vcf.gz`), number of threads.

**Outputs:** Bgzipped dosage matrix (`<vcf_basename>.dose.tsv.gz`) and its tabix index.
