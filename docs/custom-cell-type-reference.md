# Custom cell-type references

CellTypeDeconvolution and PrepareCellTypeEqtlWorkflow accept an optional File input `cell_type_mapping`. Supply a TSV with exactly two columns: `source` and `target`. Sources must uniquely and completely cover signature columns (or precomputed proportion columns). Targets identify output groups. Multiple sources may map to one target; their estimated proportions are summed after HSPE, without averaging reference profiles. Missing sources, duplicate sources, blank labels, and unresolved cloud paths fail validation.

If the input is absent, the original 22-column LM22 validation and grouping remain unchanged. The signature input retains its existing name `lm22` for compatibility. Its first column must be `gene_symbol` or `Gene symbol`. Correct website-modified gene symbols before upload using an explicit, unambiguous map. Do not replace all dots with hyphens. Gene-symbol repair is separate from the cell-type grouping file.

Suggested Tabula mapping:

| source | target |
| --- | --- |
| B | B cells |
| CD4_T | CD4 T cells |
| CD8_T | CD8 T cells |
| Erythroid | Erythroid |
| Monocyte_macrophage | Monocyte/myeloid |
| Neutrophil | Neutrophils |
| NK | NK cells |
| Plasma | Plasma cells |
| Platelet | Platelets |

The combined monocyte/macrophage label remains broad. Existing external-reference comparison against monocytes carries its existing caveat. Erythroid, platelet, and separate plasma outputs have no matching lineage in the current comparison reference and receive `no_reference_cell_type`; they are not relabeled as leukocytes. BED slugs derive from target names, e.g. CD4 T cells becomes cd4_t_cells. Existing cohorts retain their established naming.

This change affects the batched HSPE workflow and proportion processing. The separate legacy run_hspe.R CLI still uses LM22 defaults. Existing hspe_marker_fraction settings still apply: supplying a website signature does not automatically reproduce the pilot's use of all genes.

Deployment: new scripts require a rebuilt estimation image and an explicit updated estimation_docker_image. Existing pinned images do not include this feature. No image was built locally and no Terra job was submitted. WDL syntax and local tests do not establish Terra validation.
