import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(os.environ["RNASEQC_MERGE_SCRIPT"])
    if "RNASEQC_MERGE_SCRIPT" in os.environ
    else Path(__file__).resolve().parents[2]
    / "scripts" / "expression" / "merge_rnaseqc.py"
)


def write_text(path: Path, text: str) -> None:
    if path.suffix == ".gz":
        with gzip.open(path, "wt", newline="") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", newline="") as handle:
            return handle.read()
    return path.read_text()


class MergeRnaseqcTest(unittest.TestCase):
    def test_combined_manifest_validation_reports_batches_without_insert_sizes(
        self
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            prefix_file = temp / "prefix.txt"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
                "GTEX-001\tgs://bucket/a/shared.gct.gz\tgs://bucket/a/shared.gct.gz\t"
                "gs://bucket/a/shared.gct.gz\tgs://bucket/a/shared.tsv\n"
                "GTEX-002\tgs://bucket/b/shared.gct.gz\tgs://bucket/b/shared.gct.gz\t"
                "gs://bucket/b/shared.gct.gz\tgs://bucket/b/shared.tsv\n"
            )
            prefix_file.write_text("cohort.v1\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "validate-manifest",
                    "--input",
                    str(manifest),
                    "--batch-size",
                    "1",
                    "--prefix-file",
                    str(prefix_file),
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((temp / "sample_count.txt").read_text(), "2\n")
            self.assertEqual((temp / "batch_count.txt").read_text(), "2\n")
            self.assertEqual((temp / "include_insert_sizes.txt").read_text(), "false\n")

    def test_combined_manifest_validation_detects_insert_sizes(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\tinsert_size_hist\n"
                "GTEX-001\tgs://bucket/a/tpm.gct\tgs://bucket/a/count.gct\t"
                "gs://bucket/a/exon.gct\tgs://bucket/a/metrics.tsv\t"
                "gs://bucket/a/fragmentSizes.txt\n"
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "validate-manifest",
                    "--input",
                    str(manifest),
                    "--batch-size",
                    "100",
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((temp / "include_insert_sizes.txt").read_text(), "true\n")

    def test_combined_manifest_rejects_unsafe_sample_ids_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
                "sample/one\tgs://bucket/tpm.gct\tgs://bucket/count.gct\t"
                "gs://bucket/exon.gct\tgs://bucket/metrics.tsv\n"
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "validate-manifest",
                    "--input",
                    str(manifest),
                    "--batch-size",
                    "100",
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("invalid sample ID", result.stderr)
            self.assertFalse((temp / "sample_count.txt").exists())
            self.assertFalse((temp / "batch_count.txt").exists())
            self.assertFalse((temp / "include_insert_sizes.txt").exists())

    def test_combined_manifest_rejects_gcs_uri_without_an_object(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
                "GTEX-001\tgs://bucket\tgs://bucket/count.gct\t"
                "gs://bucket/exon.gct\tgs://bucket/metrics.tsv\n"
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "validate-manifest",
                    "--input",
                    str(manifest),
                    "--batch-size",
                    "100",
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("valid gs://bucket/object URI", result.stderr)
            self.assertFalse((temp / "sample_count.txt").exists())

    def test_combined_manifest_rejects_unsafe_output_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            prefix_file = temp / "prefix.txt"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
                "GTEX-001\tgs://bucket/tpm.gct\tgs://bucket/count.gct\t"
                "gs://bucket/exon.gct\tgs://bucket/metrics.tsv\n"
            )
            prefix_file.write_text("bad$(touch pwned)\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "validate-manifest",
                    "--input",
                    str(manifest),
                    "--batch-size",
                    "100",
                    "--prefix-file",
                    str(prefix_file),
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("invalid output prefix", result.stderr)
            self.assertFalse((temp / "pwned").exists())
            self.assertFalse((temp / "sample_count.txt").exists())

    def test_combined_manifest_rejects_hidden_output_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            prefix_file = temp / "prefix.txt"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
                "GTEX-001\tgs://bucket/tpm.gct\tgs://bucket/count.gct\t"
                "gs://bucket/exon.gct\tgs://bucket/metrics.tsv\n"
            )
            prefix_file.write_text(".hidden\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "validate-manifest",
                    "--input",
                    str(manifest),
                    "--batch-size",
                    "100",
                    "--prefix-file",
                    str(prefix_file),
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("invalid output prefix", result.stderr)
            self.assertFalse((temp / "sample_count.txt").exists())

    def test_prepare_batch_uses_sample_specific_names_for_shared_basenames(self) -> (
        None
    ):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            manifest.write_text(
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
                "GTEX-001\tgs://bucket/a/shared.gct.gz\tgs://bucket/a/shared.gct.gz\t"
                "gs://bucket/a/shared.gct.gz\tgs://bucket/a/shared.tsv\n"
                "GTEX-002\tgs://bucket/b/shared.gct.gz\tgs://bucket/b/shared.gct.gz\t"
                "gs://bucket/b/shared.gct.gz\tgs://bucket/b/shared.tsv\n"
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "prepare-batch",
                    "--input",
                    str(manifest),
                    "--batch-index",
                    "0",
                    "--batch-size",
                    "100",
                    "--staging-directory",
                    "individual_outputs",
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (temp / "batch_sample_ids.list").read_text(), "GTEX-001\nGTEX-002\n"
            )
            tpm_destinations = (temp / "local_tpm.list").read_text().splitlines()
            self.assertEqual(
                tpm_destinations,
                [
                    "individual_outputs/000001.GTEX-001.gene_tpm.gct.gz",
                    "individual_outputs/000002.GTEX-002.gene_tpm.gct.gz",
                ],
            )
            transfer_rows = (temp / "transfers.tsv").read_text().splitlines()
            self.assertIn(
                "gs://bucket/a/shared.gct.gz\t"
                "individual_outputs/000001.GTEX-001.gene_tpm.gct.gz",
                transfer_rows,
            )
            self.assertIn(
                "gs://bucket/b/shared.gct.gz\t"
                "individual_outputs/000002.GTEX-002.gene_tpm.gct.gz",
                transfer_rows,
            )

    def test_prepare_batch_can_skip_exon_transfers_and_keep_qc(self) -> None:
        for inserts in (False, True):
            with (self.subTest(inserts=inserts),
                  tempfile.TemporaryDirectory() as temp_dir):
                temp = Path(temp_dir)
                manifest = temp / "samples.tsv"
                header = "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv"
                row = (
                    "sample_1\tgs://test/tpm.gct\tgs://test/count.gct\t"
                    "gs://test/unused-exon.gct\tgs://test/metrics.tsv"
                )
                if inserts:
                    header += "\tinsert_size_hist"
                    row += "\tgs://test/insert.txt"
                manifest.write_text(header + "\n" + row + "\n")
                result = subprocess.run(
                    [sys.executable, str(SCRIPT), "prepare-batch", "--input",
                     str(manifest), "--batch-index", "0", "--batch-size", "100",
                     "--staging-directory", "individual_outputs", "--skip-exons"],
                    cwd=temp, capture_output=True, text=True, check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                transfers = (temp / "transfers.tsv").read_text().splitlines()
                sources = [line.split("\t")[0] for line in transfers]
                expected = ["gs://test/tpm.gct", "gs://test/count.gct",
                            "gs://test/metrics.tsv"]
                if inserts:
                    expected.append("gs://test/insert.txt")
                self.assertEqual(sources, expected)
                self.assertFalse((temp / "local_exon.list").exists())
                self.assertEqual(
                    (temp / "local_metrics.list").read_text(),
                    "individual_outputs/000001.sample_1.metrics.tsv\n",
                )

    def test_prepare_later_batch_preserves_order_and_insert_size_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            manifest = temp / "samples.tsv"
            header = (
                "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\t"
                "insert_size_hist\n"
            )
            rows = []
            for sample_number in range(1, 4):
                sample_id = f"GTEX-{sample_number:03d}"
                root = f"gs://bucket/{sample_id}"
                rows.append(
                    f"{sample_id}\t{root}/tpm.gct.gz\t{root}/count.gct.gz\t"
                    f"{root}/exon.gct.gz\t{root}/metrics.tsv.gz\t"
                    f"{root}/fragmentSizes.txt.gz\n"
                )
            manifest.write_text(header + "".join(rows))

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "prepare-batch",
                    "--input",
                    str(manifest),
                    "--batch-index",
                    "1",
                    "--batch-size",
                    "2",
                    "--staging-directory",
                    "individual_outputs",
                ],
                cwd=temp,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (temp / "batch_sample_ids.list").read_text(), "GTEX-003\n"
            )
            self.assertEqual(
                (temp / "local_insert.list").read_text(),
                "individual_outputs/000003.GTEX-003.fragmentSizes.txt.gz\n",
            )

    def test_gct_merge_streams_every_sample_column(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch_1 = temp / "batch_1.gct.gz"
            batch_2 = temp / "batch_2.gct"
            input_list = temp / "gcts.list"
            output = temp / "merged.gct.gz"
            samples = temp / "samples.txt"

            write_text(
                batch_1,
                "#1.2\n"
                "2\t2\n"
                "Name\tDescription\tsample_1\tsample_2\n"
                "gene_1\tGene one\t1\t2\n"
                "gene_2\tGene two\t3\t4\n",
            )
            write_text(
                batch_2,
                "#1.2\n"
                "2\t1\n"
                "Name\tDescription\tsample_3\n"
                "gene_1\tGene one\t5\n"
                "gene_2\tGene two\t6\n",
            )
            input_list.write_text(f"{batch_1}\n{batch_2}\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--output",
                    str(output),
                    "--sample-output",
                    str(samples),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "#1.2\n"
                "2\t3\n"
                "Name\tDescription\tsample_1\tsample_2\tsample_3\n"
                "gene_1\tGene one\t1\t2\t5\n"
                "gene_2\tGene two\t3\t4\t6\n",
            )
            self.assertEqual(samples.read_text(), "sample_1\nsample_2\nsample_3\n")

    def test_two_level_gct_merge_equals_one_level_merge(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            individual_paths = []
            for sample_number in range(1, 5):
                path = temp / f"sample_{sample_number}.gct"
                write_text(
                    path,
                    "#1.2\n"
                    "2\t1\n"
                    f"Name\tDescription\tsample_{sample_number}\n"
                    f"gene_1\tGene one\t{sample_number}\n"
                    f"gene_2\tGene two\t{sample_number * 10}\n",
                )
                individual_paths.append(path)

            direct_list = temp / "direct.list"
            direct_list.write_text("".join(f"{path}\n" for path in individual_paths))
            direct_output = temp / "direct.gct.gz"
            direct_result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(direct_list),
                    "--output",
                    str(direct_output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(direct_result.returncode, 0, direct_result.stderr)

            batch_paths = []
            for batch_number, batch_members in enumerate(
                (individual_paths[:2], individual_paths[2:]), 1
            ):
                batch_list = temp / f"batch_{batch_number}.list"
                batch_list.write_text("".join(f"{path}\n" for path in batch_members))
                batch_output = temp / f"batch_{batch_number}.gct.gz"
                batch_result = subprocess.run(
                    [
                        sys.executable,
                        str(SCRIPT),
                        "gct",
                        "--input-list",
                        str(batch_list),
                        "--output",
                        str(batch_output),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(batch_result.returncode, 0, batch_result.stderr)
                batch_paths.append(batch_output)

            batch_list = temp / "batches.list"
            batch_list.write_text("".join(f"{path}\n" for path in batch_paths))
            hierarchical_output = temp / "hierarchical.gct.gz"
            hierarchical_result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(batch_list),
                    "--output",
                    str(hierarchical_output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(
                hierarchical_result.returncode, 0, hierarchical_result.stderr
            )
            self.assertEqual(read_text(hierarchical_output), read_text(direct_output))

    def test_final_gct_merge_preserves_9042_samples_from_91_batches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch_paths = []
            next_sample_number = 1

            for batch_number in range(1, 92):
                batch_size = 100 if batch_number <= 90 else 42
                sample_names = [
                    f"sample_{sample_number:05d}"
                    for sample_number in range(
                        next_sample_number, next_sample_number + batch_size
                    )
                ]
                batch_path = temp / f"batch_{batch_number:03d}.gct.gz"
                write_text(
                    batch_path,
                    "#1.2\n" f"1\t{batch_size}\n"
                    + "\t".join(("Name", "Description", *sample_names))
                    + "\n"
                    + "\t".join(("gene_1", "Gene one", *("1" for _ in sample_names)))
                    + "\n",
                )
                batch_paths.append(batch_path)
                next_sample_number += batch_size

            input_list = temp / "batches.list"
            input_list.write_text("".join(f"{path}\n" for path in batch_paths))
            output = temp / "cohort.gct.gz"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            output_lines = read_text(output).splitlines()
            self.assertEqual(output_lines[1], "1\t9042")
            self.assertEqual(len(output_lines[2].split("\t")), 9044)
            self.assertEqual(len(output_lines[3].split("\t")), 9044)
            self.assertTrue(output_lines[2].endswith("\tsample_09042"))

    def test_gct_merge_rejects_a_late_feature_mismatch_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch_1 = temp / "batch_1.gct"
            batch_2 = temp / "batch_2.gct"
            input_list = temp / "batches.list"
            output = temp / "cohort.gct.gz"

            common_prefix = "#1.2\n2\t1\nName\tDescription\t{}\ngene_1\tGene one\t1\n"
            write_text(
                batch_1,
                common_prefix.format("sample_1") + "gene_2\tGene two\t2\n",
            )
            write_text(
                batch_2,
                common_prefix.format("sample_2") + "wrong_gene\tGene two\t3\n",
            )
            input_list.write_text(f"{batch_1}\n{batch_2}\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("feature differs at GCT row 2", result.stderr)
            self.assertFalse(output.exists())

    def test_gct_merge_rejects_duplicate_sample_names(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch_1 = temp / "batch_1.gct"
            batch_2 = temp / "batch_2.gct"
            input_list = temp / "batches.list"
            output = temp / "cohort.gct.gz"

            gct_text = (
                "#1.2\n"
                "1\t1\n"
                "Name\tDescription\tsample_1\n"
                "gene_1\tGene one\t1\n"
            )
            write_text(batch_1, gct_text)
            write_text(batch_2, gct_text)
            input_list.write_text(f"{batch_1}\n{batch_2}\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("Duplicate sample names", result.stderr)
            self.assertFalse(output.exists())

    def test_gct_merge_rejects_an_unexpected_sample_order(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            gct = temp / "counts.gct"
            input_list = temp / "gcts.list"
            expected_samples = temp / "samples.txt"
            output = temp / "merged.gct.gz"

            write_text(
                gct,
                "#1.2\n"
                "1\t2\n"
                "Name\tDescription\tsample_2\tsample_1\n"
                "gene_1\tGene one\t2\t1\n",
            )
            input_list.write_text(f"{gct}\n")
            expected_samples.write_text("sample_1\nsample_2\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--expected-samples",
                    str(expected_samples),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("sample names do not match", result.stderr)
            self.assertFalse(output.exists())

    def test_gct_merge_can_assign_sample_names_from_manifest_order(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            gct_1 = temp / "000001.gct"
            gct_2 = temp / "000002.gct"
            input_list = temp / "gcts.list"
            sample_names = temp / "sample_ids.txt"
            output = temp / "merged.gct.gz"

            common_header = "#1.2\n1\t1\nName\tDescription\ttransqtl\n"
            write_text(gct_1, common_header + "gene_1\tGene one\t1\n")
            write_text(gct_2, common_header + "gene_1\tGene one\t2\n")
            input_list.write_text(f"{gct_1}\n{gct_2}\n")
            sample_names.write_text("sample_1\nsample_2\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--sample-names",
                    str(sample_names),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "#1.2\n"
                "1\t2\n"
                "Name\tDescription\tsample_1\tsample_2\n"
                "gene_1\tGene one\t1\t2\n",
            )

    def test_gct_merge_rejects_tab_in_assigned_sample_name(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            gct = temp / "000001.gct"
            input_list = temp / "gcts.list"
            sample_names = temp / "sample_ids.txt"
            output = temp / "merged.gct.gz"

            write_text(
                gct,
                "#1.2\n"
                "1\t1\n"
                "Name\tDescription\ttransqtl\n"
                "gene_1\tGene one\t1\n",
            )
            input_list.write_text(f"{gct}\n")
            sample_names.write_text("sample\t1\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--sample-names",
                    str(sample_names),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("invalid sample ID", result.stderr)
            self.assertFalse(output.exists())

    def test_gct_merge_rejects_surrounding_sample_name_whitespace(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            gct = temp / "000001.gct"
            input_list = temp / "gcts.list"
            sample_names = temp / "sample_ids.txt"
            output = temp / "merged.gct.gz"

            write_text(
                gct,
                "#1.2\n"
                "1\t1\n"
                "Name\tDescription\ttransqtl\n"
                "gene_1\tGene one\t1\n",
            )
            input_list.write_text(f"{gct}\n")
            sample_names.write_text(" sample_1\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "gct",
                    "--input-list",
                    str(input_list),
                    "--sample-names",
                    str(sample_names),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("invalid sample ID", result.stderr)
            self.assertFalse(output.exists())

    def test_individual_metrics_are_appended_as_sample_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            sample_1 = temp / "000001.metrics.tsv"
            sample_2 = temp / "000002.metrics.tsv"
            input_list = temp / "metrics.list"
            expected_samples = temp / "samples.txt"
            output = temp / "cohort.metrics.txt.gz"

            sample_1.write_text(
                "Sample\tsample_1\n" "Mapping Rate\t0.9\n" "Total Reads\t100\n"
            )
            sample_2.write_text(
                "Sample\tsample_2\n" "Mapping Rate\t0.8\n" "Total Reads\t200\n"
            )
            input_list.write_text(f"{sample_1}\n{sample_2}\n")
            expected_samples.write_text("sample_1\nsample_2\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "metrics-individual",
                    "--input-list",
                    str(input_list),
                    "--expected-samples",
                    str(expected_samples),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "sample_id\tMapping Rate\tTotal Reads\n"
                "sample_1\t0.9\t100\n"
                "sample_2\t0.8\t200\n",
            )

    def test_individual_metrics_can_use_assigned_sample_names(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            metrics_1 = temp / "000001.metrics.tsv"
            metrics_2 = temp / "000002.metrics.tsv"
            input_list = temp / "metrics.list"
            sample_names = temp / "sample_ids.txt"
            output = temp / "cohort.metrics.txt.gz"

            metrics_1.write_text("Sample\ttransqtl\nMapping Rate\t0.9\n")
            metrics_2.write_text("Sample\ttransqtl\nMapping Rate\t0.8\n")
            input_list.write_text(f"{metrics_1}\n{metrics_2}\n")
            sample_names.write_text("sample_1\nsample_2\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "metrics-individual",
                    "--input-list",
                    str(input_list),
                    "--sample-names",
                    str(sample_names),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "sample_id\tMapping Rate\n" "sample_1\t0.9\n" "sample_2\t0.8\n",
            )

    def test_aggregated_metrics_are_appended_without_repeating_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch_1 = temp / "batch_1.metrics.txt.gz"
            batch_2 = temp / "batch_2.metrics.txt.gz"
            input_list = temp / "metrics.list"
            expected_samples = temp / "samples.txt"
            output = temp / "cohort.metrics.txt.gz"

            write_text(
                batch_1,
                "sample_id\tMapping Rate\tTotal Reads\n"
                "sample_1\t0.9\t100\n"
                "sample_2\t0.8\t200\n",
            )
            write_text(
                batch_2,
                "sample_id\tMapping Rate\tTotal Reads\n" "sample_3\t0.7\t300\n",
            )
            input_list.write_text(f"{batch_1}\n{batch_2}\n")
            expected_samples.write_text("sample_1\nsample_2\nsample_3\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "metrics-aggregated",
                    "--input-list",
                    str(input_list),
                    "--expected-samples",
                    str(expected_samples),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "sample_id\tMapping Rate\tTotal Reads\n"
                "sample_1\t0.9\t100\n"
                "sample_2\t0.8\t200\n"
                "sample_3\t0.7\t300\n",
            )

    def test_failed_metrics_validation_does_not_publish_an_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch = temp / "batch.metrics.txt.gz"
            input_list = temp / "metrics.list"
            expected_samples = temp / "samples.txt"
            output = temp / "cohort.metrics.txt.gz"

            write_text(
                batch,
                "sample_id\tMapping Rate\n" "wrong_sample\t0.9\n",
            )
            input_list.write_text(f"{batch}\n")
            expected_samples.write_text("sample_1\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "metrics-aggregated",
                    "--input-list",
                    str(input_list),
                    "--expected-samples",
                    str(expected_samples),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("sample names do not match", result.stderr)
            self.assertFalse(output.exists())

    def test_individual_insert_sizes_use_the_union_of_bins(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            sample_1 = temp / "000001.fragmentSizes.txt"
            sample_2 = temp / "000002.fragmentSizes.txt"
            input_list = temp / "insert_sizes.list"
            expected_samples = temp / "samples.txt"
            output = temp / "cohort.insert_size_hists.txt.gz"

            sample_1.write_text("100\t5\n102\t2\n")
            sample_2.write_text("101\t7\n102\t3\n")
            input_list.write_text(f"{sample_1}\n{sample_2}\n")
            expected_samples.write_text("sample_1\nsample_2\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "insert-sizes-individual",
                    "--input-list",
                    str(input_list),
                    "--expected-samples",
                    str(expected_samples),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "insert_size\tsample_1\tsample_2\n"
                "100\t5\t0\n"
                "101\t0\t7\n"
                "102\t2\t3\n",
            )

    def test_aggregated_insert_sizes_use_a_streaming_outer_join(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            batch_1 = temp / "batch_1.insert_size_hists.txt.gz"
            batch_2 = temp / "batch_2.insert_size_hists.txt.gz"
            input_list = temp / "insert_sizes.list"
            expected_samples = temp / "samples.txt"
            output = temp / "cohort.insert_size_hists.txt.gz"

            write_text(
                batch_1,
                "insert_size\tsample_1\tsample_2\n" "100\t5\t0\n" "102\t2\t3\n",
            )
            write_text(
                batch_2,
                "insert_size\tsample_3\n" "101\t7\n" "102\t8\n",
            )
            input_list.write_text(f"{batch_1}\n{batch_2}\n")
            expected_samples.write_text("sample_1\nsample_2\nsample_3\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "insert-sizes-aggregated",
                    "--input-list",
                    str(input_list),
                    "--expected-samples",
                    str(expected_samples),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                read_text(output),
                "insert_size\tsample_1\tsample_2\tsample_3\n"
                "100\t5\t0\t0\n"
                "101\t0\t0\t7\n"
                "102\t2\t3\t8\n",
            )


if __name__ == "__main__":
    unittest.main()
