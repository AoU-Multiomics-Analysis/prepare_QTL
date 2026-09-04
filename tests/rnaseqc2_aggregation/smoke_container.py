"""Run real container and WDL smoke checks on a GitHub Actions Docker host."""

import argparse
import gzip
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WDL = ROOT / "workflows/expression/rnaseqc2_aggregate_batched.wdl"
SCRIPT = "/opt/prepare_qtl/scripts/expression/merge_rnaseqc.py"


def run_task(image: str, task: str, inputs: dict, directory: Path) -> dict:
    directory.mkdir()
    inputs.update(docker_image=image, memory_gb=1, disk_space_gb=1, num_preempt=0)
    input_json = directory / "inputs.json"
    output_json = directory / "outputs.json"
    input_json.write_text(json.dumps(inputs))
    subprocess.run(
        [
            "miniwdl", "run", str(WDL), "--task", task,
            "--input", str(input_json), "--dir", str(directory / "run"),
            "-o", str(output_json), "--no-cache", "--verbose",
        ],
        check=True,
    )
    return {
        key.removeprefix(task + "."): value
        for key, value in json.loads(output_json.read_text())["outputs"].items()
    }


def read_gzip(path: str) -> str:
    with gzip.open(path, "rt") as handle:
        return handle.read()


def check_wdl_tasks(image: str, temp: Path, include_inserts: bool) -> None:
    temp.mkdir()
    manifest = temp / "samples.tsv"
    header = "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv"
    if include_inserts:
        header += "\tinsert_size_hist"
    rows = []
    for sample in ("sample_1", "sample_2", "sample_3"):
        row = [sample] + [f"gs://smoke/{sample}/shared.{suffix}" for suffix in
                          ("tpm.gct", "count.gct", "exon.gct", "metrics.tsv")]
        if include_inserts:
            row.append(f"gs://smoke/{sample}/shared.hist.txt")
        rows.append("\t".join(row))
    manifest.write_text(header + "\n" + "\n".join(rows) + "\n")
    validated = run_task(
        image, "validate_rnaseqc_manifests",
        dict(sample_manifest=str(manifest), prefix="smoke", batch_size=2),
        temp / "validate",
    )
    prefix_file = validated.pop("prefix_file")
    assert Path(prefix_file).read_text() == "smoke\n"
    assert validated == {
        "sample_count": 3, "batch_count": 2, "include_insert_sizes": include_inserts,
    }, validated

    inputs = {"prefix_file": prefix_file, "include_insert_sizes": include_inserts}
    expected_gcts = {}
    for kind, label in (("tpm", "gene_tpm"), ("count", "gene_reads"),
                        ("exon_count", "exon_reads")):
        paths = []
        for batch, columns, values in (
            (1, "sample_1\tsample_2", "1\t2"), (2, "sample_3", "3"),
        ):
            path = temp / f"batch_{batch}.{label}.gct"
            n = 2 if batch == 1 else 1
            path.write_text(
                f"#1.2\n1\t{n}\nName\tDescription\t{columns}\n"
                f"feature_1\tFeature one\t{values}\n"
            )
            paths.append(str(path))
        inputs[f"batch_{kind}_gcts"] = paths
        expected_gcts[f"{kind}_gct"] = (
            "#1.2\n1\t3\nName\tDescription\tsample_1\tsample_2\tsample_3\n"
            "feature_1\tFeature one\t1\t2\t3\n"
        )
    metrics = []
    for batch, rows in ((1, "sample_1\t10\nsample_2\t20\n"),
                        (2, "sample_3\t30\n")):
        path = temp / f"batch_{batch}.metrics.txt"
        path.write_text("sample_id\tTotal Reads\n" + rows)
        metrics.append(str(path))
    inputs["batch_metrics"] = metrics
    inputs["batch_insert_size_hists"] = []
    if include_inserts:
        for batch, text in (
            (1, "insert_size\tsample_1\tsample_2\n100\t1\t2\n"),
            (2, "insert_size\tsample_3\n101\t3\n"),
        ):
            path = temp / f"batch_{batch}.insert_size_hists.txt"
            path.write_text(text)
            inputs["batch_insert_size_hists"].append(str(path))
    outputs = run_task(image, "merge_rnaseqc_batches", inputs, temp / "merge")
    for name, expected in expected_gcts.items():
        assert read_gzip(outputs[name]) == expected, name
    assert read_gzip(outputs["metrics"]) == (
        "sample_id\tTotal Reads\nsample_1\t10\nsample_2\t20\nsample_3\t30\n"
    )
    if include_inserts:
        assert len(outputs["insert_size_hists"]) == 1
        assert read_gzip(outputs["insert_size_hists"][0]) == (
            "insert_size\tsample_1\tsample_2\tsample_3\n"
            "100\t1\t2\t0\n101\t0\t0\t3\n"
        )
    else:
        assert outputs["insert_size_hists"] == []


def check_workflow_rejects_unsafe_prefix(
    image: str, temp: Path, manifest: Path,
) -> None:
    """Start the whole workflow; validation must fail before any GCS transfer."""
    temp.mkdir()
    input_json = temp / "inputs.json"
    workflow = "rnaseqc2_aggregate_batched_workflow."
    input_json.write_text(json.dumps({
        workflow + "sample_manifest": str(manifest),
        workflow + "prefix": "x$(touch prefix_injection_marker)",
        workflow + "docker_image": image,
        workflow + "merge_disk_space_gb": 1,
        workflow + "num_preempt": 0,
    }))
    result = subprocess.run(
        [
            "miniwdl", "run", str(WDL), "--input", str(input_json),
            "--dir", str(temp / "run"), "--no-cache", "--verbose",
        ],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, "The workflow accepted an unsafe prefix"
    task_errors = "\n".join(path.read_text() for path in temp.rglob("stderr.txt"))
    assert "invalid output prefix" in task_errors, result.stdout + result.stderr
    assert not list(temp.rglob("prefix_injection_marker")), "Prefix ran as shell code"
    assert not list(temp.rglob("prefix.txt")), "Invalid prefix was published"
    print("The workflow rejected an unsafe prefix in the validation task.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    args = parser.parse_args()
    # Mount only tests. The production script must come from the image.
    # Bypass the entrypoint, as a WDL backend can do, to test the runtime PATH.
    subprocess.run(
        [
            "docker", "run", "--rm", "--network", "none",
            "--entrypoint", "/bin/bash",
            "--volume", f"{ROOT / 'tests/rnaseqc2_aggregation'}:/tests:ro",
            "--workdir", "/tmp", "--env", f"RNASEQC_MERGE_SCRIPT={SCRIPT}",
            args.image, "-c", r"""
set -euo pipefail
python3 -m unittest discover -s /tests -p test_merge_rnaseqc.py -v
python3 /tests/check_gce_credentials.py
python3 - <<'PY'
import subprocess
text = subprocess.check_output(["gsutil", "version", "-l"], text=True)
print(text)
assert "compiled crcmod: True" in text
PY
printf 'copy-check\n' > source.txt
printf '%s\0%s\0' source.txt copied.txt |
    xargs -0 -n 2 -P 2 bash -c 'gsutil cp "$1" "$2"' _
test "$(awk 'END {print NR}' copied.txt)" = 1
python3 - <<'PY'
from pathlib import Path
assert Path("copied.txt").read_text() == "copy-check\n"
PY
date -u '+%Y-%m-%dT%H:%M:%SZ'
""",
        ],
        check=True,
    )
    with tempfile.TemporaryDirectory(prefix="rnaseqc-container-") as temp_dir:
        temp = Path(temp_dir)
        for include_inserts in (False, True):
            check_wdl_tasks(args.image, temp / str(include_inserts), include_inserts)
        check_workflow_rejects_unsafe_prefix(
            args.image, temp / "unsafe_prefix", temp / "False" / "samples.tsv",
        )
    print("Container and WDL smoke checks passed.")


if __name__ == "__main__":
    main()
