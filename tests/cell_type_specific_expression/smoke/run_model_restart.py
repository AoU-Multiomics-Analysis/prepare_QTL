"""Reuse a smoke model and verify export/QTL outputs without a fit."""
import argparse
import gzip
import json
from pathlib import Path
import subprocess


PREFIX = "PrepareCellTypeEqtlWorkflow."
WORKFLOW = "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--deconvolution-only', action='store_true',
                        help='Check model reuse and BED export without repeating QTL preparation')
    parser.add_argument(
        "--baseline-outputs",
        type=Path,
        default=Path("ci-runs/precomputed/outputs.json"),
        help="Flat MiniWDL outputs.json from the completed baseline run",
    )
    parser.add_argument(
        "--baseline-inputs",
        type=Path,
        default=Path(
            "tests/cell_type_specific_expression/fixtures/precomputed-e2e.inputs.json"
        ),
        help="Input JSON used for the completed baseline run",
    )
    parser.add_argument(
        "--restart-inputs",
        type=Path,
        default=Path("ci-runs/model-restart.inputs.json"),
        help="Path for the generated restart input JSON",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("ci-runs/model-restart"),
        help="Exact MiniWDL output directory for the restart run",
    )
    return parser.parse_args(argv)


def read_outputs(path):
    # MiniWDL writes a flat object. Its CLI stdout alone uses an outputs wrapper.
    return json.loads(path.read_text())


def main(argv=None):
    args = parse_args(argv)
    baseline = read_outputs(args.baseline_outputs)
    inputs = json.loads(args.baseline_inputs.read_text())
    inputs[PREFIX + "precomputed_tca_model"] = baseline[PREFIX + "tca_model"]
    # Both files exist, but are invalid as proportions/covariates. They must be ignored.
    inputs[PREFIX + "precomputed_proportions"] = inputs[PREFIX + "lm22"]
    inputs[PREFIX + "deconvolution_covariates"] = inputs[PREFIX + "lm22"]
    inputs[PREFIX + "gene_type"] = ["ignored_for_model_restart"]
    inputs[PREFIX + "OutputPrefix"] = "synthetic.restart"
    workflow = WORKFLOW
    if args.deconvolution_only:
        import WDL
        workflow = 'workflows/cell_type_specific_expression/deconvolution.wdl'
        allowed = {d.name for d in WDL.load(workflow).workflow.inputs}
        names = {key.removeprefix(PREFIX): value for key, value in inputs.items()}
        if 'deconvolution_covariates' in names:
            names['covariates'] = names.pop('deconvolution_covariates')
        inputs = {'CellTypeDeconvolution.' + key: value for key, value in names.items()
                  if key in allowed}
    args.restart_inputs.parent.mkdir(parents=True, exist_ok=True)
    args.restart_inputs.write_text(json.dumps(inputs, indent=2))
    subprocess.run(
        [
            "miniwdl",
            "run",
            workflow,
            "--input",
            str(args.restart_inputs),
            "--dir",
            str(args.output_directory) + "/.",
            "--verbose",
            "--no-color",
        ],
        check=True,
    )
    result = read_outputs(args.output_directory / "outputs.json")
    if args.deconvolution_only:
        result = {PREFIX + key.removeprefix('CellTypeDeconvolution.'): value
                  for key, value in result.items()}
    for name in (
        "estimated_proportions",
        "tca_model_unfiltered",
        "fit_tca_log",
        "proportions_lm22",
        "proportions_combined",
        "gene_type_filter_log",
    ):
        assert result[PREFIX + name] is None, f"{name} should be skipped"
    for name in ("cell_type_beds", "filtered_cell_type_beds"):
        expected = baseline[PREFIX + name]
        actual = result[PREFIX + name]
        assert len(actual) == len(expected) > 0
        for left, right in zip(expected, actual):
            with gzip.open(left, "rt") as baseline_bed, gzip.open(
                right, "rt"
            ) as restart_bed:
                assert (
                    baseline_bed.read() == restart_bed.read()
                ), f"Restarted {name} differs from the original export"
    qtl_outputs = (
        "int_beds",
        "scaled_beds",
        "int_phenotype_pcs",
        "scaled_phenotype_pcs",
        "int_merged_covariates",
        "scaled_merged_covariates",
    )
    for name in (() if args.deconvolution_only else qtl_outputs):
        assert len(result[PREFIX + name]) == len(result[PREFIX + "cell_type_beds"])
        assert all(Path(path).is_file() for path in result[PREFIX + name])
    if not args.deconvolution_only:
        assert Path(result[PREFIX + "cell_type_qtl_manifest"]).is_file()
    expected_images = dict(baseline[PREFIX + 'stage_images'])
    if args.deconvolution_only:
        expected_images.pop('qtl')
    assert result[PREFIX + "stage_images"] == expected_images
    parameters = json.loads(
        Path(result[PREFIX + "effective_parameters_file"]).read_text()
    )
    assert parameters["proportion_mode"] == "precomputed_model"
    assert (
        "tca_max_iters" not in parameters
    ), "Do not report new fitting settings for an old model"
    print(
        "Model restart: skipped fit; identical exported BEDs; "
        + ("deconvolution outputs present." if args.deconvolution_only else "QTL outputs and manifest present.")
    )


if __name__ == "__main__":
    main()
