# Trans-QTL regions from ancestry-specific LD

This WDL 1.0 workflow defines conservative candidate fine-mapping regions for Terra. It uses selected trans associations, PLINK genotypes, and sample ancestry assignments. It does not run fine-mapping.

Every region has a minimum width of **2 Mb total**. Regions can be wider. The workflow combines the earliest and latest LD-supported boundaries across ancestry groups, adds padding, extends short intervals, and merges overlaps. All source phenotype associations remain in the outputs.

## Region rules

1. Read the association tables and retain rows with `p_value <= pvalue_threshold`. The input must already use the intended trans definition and significance procedure. This workflow does not calculate FDR, correct multiple testing, or remove cis associations.
2. Use every distinct retained associated variant as a starting variant. Do not select only one lead variant or prune away signals from other phenotypes.
3. Select samples separately for each ancestry label. Use only samples recorded as founders in FAM/PSAM. Retain every label in the assignment file, including MID and ADMIX. Founder status does not establish that samples are unrelated; supply an appropriate LD sample subset if needed.
4. Calculate allele counts in that selected sample group. Retain biallelic variants with **MAF strictly greater than 1%** by default. Apply the rule to both members of an LD pair. There is **no MAC cutoff**. A variant need not pass the frequency rule in every group.
5. Calculate PLINK 2 `--r2-unphased` between starting variants and eligible nearby genotype variants. An LD partner does not need to appear in the association table. The default reporting threshold is `r² >= 0.1`.
6. Start the LD search at 1 Mb on each side. If a qualifying partner is within the outer 100 kb of the search boundary, double the search distance for that seed. Stop at a default maximum of 10 Mb per side. Search-limit hits retain the full searched interval and get an unresolved flag.
7. For each seed, combine the **minimum start and maximum end** across assessed groups. Add 100 kb on each side. Extend any shorter interval to **2,000,000 bp**, symmetrically where possible. Near a chromosome end, shift the extension to the available side. If the entire chromosome is shorter, return the whole chromosome and flag it.
8. If no ancestry provides qualifying LD support for a seed, retain a fallback of 1 Mb on each side, plus padding. The default fallback therefore spans up to 2.2 Mb. Apply the same minimum-width rule.
9. Merge overlapping intervals on the same chromosome. Preserve each seed, source association, phenotype, dataset, cell type, and modality. A merged region can contain multiple signals.

The LD threshold, search distances, padding, fallback and minimum width are configurable heuristics. These are **LD-defined candidate regions**, not formal haplotype blocks. Expansion checks LD with the original associated variants; it does not repeatedly promote every linked genotype variant into a new seed. A disconnected long-range LD segment beyond the search radius can be missed. Set a larger `initial_search_bp` when a wider first search is needed; setting it equal to `max_search_bp` searches that full distance immediately. Separate output regions are not guaranteed to be LD-independent.

No MAC rule suppresses small ancestry groups. The output records which ancestry groups and variant pairs support the outer boundaries, their r² values, and allele counts. A broad interval supported only by MID or ADMIX is retained. Splitting by global ancestry label does not remove all ancestry variation within an admixed group.

## Terra inputs

Upload or register the single `workflows/genotype/trans_ld_regions.wdl` file. It has no WDL imports. The runtime scripts are installed in the container described below.

Use `examples/trans_ld_regions/terra.inputs.json` as a submission template. Replace all example bucket and image values. The JSON is a Terra submission file, not a script argument wrapper.

### Genotypes

`genotypes` is an array of structs. Each entry declares a chromosome and three **File** inputs:

```json
{
  "chromosome": "1",
  "format": "bed",
  "genotype": "gs://bucket/chr1.bed",
  "variants": "gs://bucket/chr1.bim",
  "samples": "gs://bucket/chr1.fam"
}
```

For PLINK 2, use `format: "pgen"` and supply `.pgen`, `.pvar` or `.pvar.zst`, and `.psam`. PGEN files must have an embedded index; external `.pgi` files are not part of this interface. Filenames and storage directories need not share a prefix. Cromwell localizes all three files separately; the command passes their local paths directly to PLINK.

Use one entry per chromosome. Supported chromosome names are `1`–`22` and `X`, with an optional `chr` prefix. Standard BIM code `23` is recognized as X. Other chromosomes are rejected for association inputs. Chromosome X uses PLINK's unphased dosage calculation and requires correct sex metadata; split pseudoautosomal-region codes are not supported in this interface. The real X regression uses female diploid genotypes; mixed-sex X data have not been tested here.

The variant IDs must match the association table exactly. Duplicate IDs within a chromosome cause failure. Positions must agree. A missing seed is retained as unresolved, but a shard containing no variants on its declared chromosome causes failure. All retained association chromosomes must have a genotype entry.

### Ancestry samples

Supply a TSV with these exact headers:

```text
IID         ancestry
sample_001  EUR
sample_002  MID
sample_003  ADMIX
```

Use tabs between columns. IID must be unique in both the assignment table and each genotype sample file. Extra table columns are ignored. Each participant belongs to one group. Participants absent from the assignment file are excluded from LD. The QC output reports assigned participants missing from each genotype shard and excluded nonfounders. Ancestry labels are not inferred.

Groups with fewer than 50 usable founders are marked `insufficient_samples`. This is a minimum sample requirement for the LD calculation, not a minor-allele-count filter. No PLINK `--bad-ld` override is used.

### Associations

`association_files` accepts explicit Files in TSV, TSV.gz or Parquet format. All tables use the same selected column mapping. Parquet is read in batches. The default columns are:

| Column | Meaning |
|---|---|
| `variant_id` | Exact genotype variant ID |
| `phenotype_id` | Target gene, splice event or other molecular trait |
| `pval` | Numeric association p-value |
| `chrom`, `pos` | Chromosome and 1-based position; provide both or neither |
| `modality` | Optional; otherwise `default_modality` is used |
| `dataset_id` | Optional; otherwise `input_1`, `input_2`, etc. identify input files |
| `cell_type` | Optional; otherwise `unspecified` |

Each name can be changed with the corresponding workflow `*_column` input. If chromosome and position columns are both absent, IDs must follow `chr1_123_A_G` or `1:123:A:G`. IDs such as rsIDs require explicit coordinates. Data from different genome builds must not be mixed. `genome_build` is recorded as metadata; it does not perform liftover or independently identify a file's build.

For already filtered associations, use `pvalue_threshold: 1.0` to retain every valid input row. Do not use that setting on a full nominal association scan unless every row is intended to define a region. No significance threshold is inferred from the source data.

### Chromosome sizes and resources

`chrom_sizes` is a two-column chromosome/length text file or a FASTA `.fai` for the declared build. Non-supported contigs in a `.fai` are ignored. Sizes prevent padding beyond chromosome ends.

The workflow scatters over chromosome entries and ancestry groups. LD tasks default to 4 CPUs, 16 GiB RAM and 100 GB local disk. Size the disk to hold the localized genotype shard and temporary LD reports. Separate preparation and merge resources are configurable. The default maximum search is a distance limit; there is no additional limit on the final merged interval width.

## Outputs

| Output | Contents |
|---|---|
| `regions.tsv` | Region ID, chromosome, start, end, width, seed count, build and flags |
| `regions.bed` | Four-column BED for region extraction |
| `region_associations.tsv.gz` | Every retained source association, assigned to a region |
| `trans_window_associations.tsv.gz` | Region–phenotype links, compatible column names for the existing preparation pipeline |
| `seed_regions.tsv` | Each seed's padded interval, final region ID and groups supporting its left/right boundaries |
| `ancestry_intervals.tsv.gz` | Per-seed, per-ancestry LD bounds, boundary partners, r² and allele counts |
| `region_qc.json` | Region counts, fallback/search-limit counts and merge parameters |
| `preparation_qc` | Association counts, ancestry assignment counts and input settings |
| `ancestry_qc` | Sample counts, PLINK version, filter settings and status counts per task |
| `ancestry_allele_counts` | Per-variant MAF, MAC, observed allele count and eligibility per task |
| Logs | Preparation, PLINK commands/output and merge logs |

`regions.tsv`, `regions.bed`, `seed_regions.tsv`, and the exported association interval columns use **0-based, half-open** coordinates. `pos` and raw `ancestry_intervals` start/end use **1-based, inclusive** coordinates, as in PLINK. BED and compatibility output chromosome names have a `chr` prefix; other tables use `1`–`22`/`X`.

The compatibility table has `window_id`, `chrom`, `start`, `end`, `modality`, `molecular_trait_id`, `p_value`, `dataset_id`, and `cell_type`. Its p-value is the minimum across associated variants for each region/dataset/cell-type/modality/phenotype combination. Filter to the appropriate dataset and cell type before passing it to a downstream pipeline that does not model those fields. The complete row-level evidence is retained separately.

Do not use `ancestry_allele_counts` or the LD partner list as the downstream variant inclusion list. The next pipeline should extract the genomic interval and apply its own fine-mapping variant QC.

### Status meanings

- `ok`: at least one qualifying finite LD pair was observed, and no search edge remains unresolved.
- `search_limit`: qualifying LD reaches the search edge at the maximum distance; full search bounds are retained.
- `missing_variant`: seed ID is absent from the supplied genotype chromosome.
- `maf_filtered`: seed does not exceed the group's MAF threshold, or has no observed alleles.
- `multiallelic`: seed is not eligible for the biallelic LD calculation.
- `insufficient_samples`: fewer than 50 usable group founders.
- `no_reported_ld`: no qualifying LD pairs were reported. This includes correlations below the threshold, no eligible partners and undefined correlations, such as constant dosage. It is **not evidence of independence**.

Only `ok` and `search_limit` contribute measured boundaries. If some groups are unassessed, the seed receives `some_groups_unassessed`. If all groups are unassessed, it receives a fallback. No seed is silently removed for these statuses.

## Runtime image and GitHub Actions

`tools/trans_ld_regions/Dockerfile` uses `mambaorg/micromamba:2.0.5` and pinned conda-forge/bioconda packages:

- Python 3.11.11
- PLINK 2.0.0-a.6.9, Bioconda package `2.0.0a.6.9=h9948957_0`
- PyArrow 19.0.1
- zstandard 0.23.0

The repository includes dedicated Actions workflows for this tool. `trans-ld-regions-check.yml` validates WDL with miniwdl and Cromwell womtool and runs unit regressions. `trans-ld-regions-image.yml` builds a Linux amd64 image, runs real PLINK smoke tests for BED/PGEN and executes all three rendered WDL commands inside the image. Only after those checks does it publish a commit-tagged image to GHCR on the default branch as `ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions:sha-<commit>`.

Automatic image builds trigger only when a contained script under `tools/trans_ld_regions/scripts/` changes. Dockerfile or environment edits require a manual `workflow_dispatch` build. No local Docker build is needed. Supply the published image URL, preferably its digest, as `docker_image`. No image has been published by this local task.

All files opened by scripts remain typed WDL Files until command rendering. File arrays become task-local newline lists at that point. Generated `write_lines` results stay File-typed through final command-placeholder localization; they are not wrapped in `sub` or string conversion. Filenames containing newlines are not supported. Required inputs are checked for readability, and unresolved cloud URIs produce a localization error before computation. JSON files are used only for structured output/QC and Terra submission inputs.

## Verification

From the repository root, with Python 3.11, miniwdl, PyArrow and zstandard available:

```bash
miniwdl check workflows/genotype/trans_ld_regions.wdl
PLINK2=/path/to/plink2 python3 -m unittest discover -s tests/trans_ld_regions -v
java -jar /path/to/womtool-87.jar validate workflows/genotype/trans_ld_regions.wdl
```

Without `PLINK2` or PLINK on PATH, tests that need real PLINK are skipped. The GitHub image tests set `PLINK2` explicitly. See `docs/trans_ld_regions/validation.md` for the results from this implementation.

**The complete workflow has not been tested on Terra.** Local rendered-command tests simulate cloud-to-local File handling; they do not establish a successful Terra run. GitHub Actions have been supplied but have not been run for this new bundle.

## Method references

- [PLINK 2 LD calculations](https://www.cog-genomics.org/plink/2.0/ld)
- [PLINK allele count and frequency filters](https://www.cog-genomics.org/plink/2.0/filter#maf)
- [PLINK output formats](https://www.cog-genomics.org/plink/2.0/formats)

The ancestry-union, padding and minimum-width rules are the study-specific design chosen for this pipeline; they are not presented as a published standard.
