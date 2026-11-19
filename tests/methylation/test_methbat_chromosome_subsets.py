#!/usr/bin/env python3
"""Check the WDL and exercise its commands with simulated Terra localization."""

import csv
import gzip
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import WDL


ROOT = Path(__file__).resolve().parents[2]
DESCRIPTOR = ROOT / "workflows/methylation/MethBat_jointseg_chromosome_subsets.wdl"


class CloudGeneratedFiles(WDL.StdLib.Base):
    def __init__(self, directory):
        super().__init__("1.0", str(directory))
        self.paths = {}

    def _virtualize_filename(self, filename):
        uri = f"gs://execution-outputs/generated/{Path(filename).name}"
        self.paths[uri] = filename
        return uri


class WorkflowTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow = WDL.load(str(DESCRIPTOR))
        cls.tasks = {task.name: task for task in cls.workflow.tasks}

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.env = dict(os.environ, PATH=f"{self.bin}:{os.environ['PATH']}")
        if not os.environ.get("METHBAT_REAL"):
            for tool in ("methbat_prepare_pileups.py", "methbat_compare_methods.py"):
                (self.bin / tool).symlink_to(ROOT / "tools/methbat" / tool)

    def run_task(self, task_name, directory, inputs, localized, want_error=None):
        directory.mkdir(exist_ok=True)
        task = self.tasks[task_name]
        values = WDL.values_from_json(inputs, task.available_inputs)
        values = WDL.Value.rewrite_env_paths(values, lambda file: localized.get(file.value, file.value))
        stdlib = CloudGeneratedFiles(directory)
        for declaration in task.postinputs:
            values = values.bind(declaration.name, declaration.expr.eval(values, stdlib))

        original_eval = WDL.Expr.Apply._eval

        def eval_with_generated_file_mapping(expression, env, library):
            value = original_eval(expression, env, library)
            # Model localization only if a writer's File reaches the placeholder.
            if (isinstance(expression.parent, WDL.Expr.Placeholder)
                    and isinstance(value, WDL.Value.File) and value.value in library.paths):
                return WDL.Value.File(library.paths[value.value])
            return value

        with patch.object(WDL.Expr.Apply, "_eval", eval_with_generated_file_mapping):
            command = task.command.eval(values, stdlib).value
        result = subprocess.run(
            ["bash", "-c", command], cwd=directory, env=self.env, capture_output=True, text=True
        )
        if want_error:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(want_error, result.stderr)
        else:
            self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def test_no_workflow_scope_file_writes(self):
        self.assertEqual(self.workflow.wdl_version, "1.0")
        pending = [self.workflow.workflow]
        while pending:
            node = pending.pop()
            if isinstance(node, WDL.Expr.Apply):
                self.assertFalse(node.function_name.startswith("write_"), node.pos)
            pending.extend(node.children)

    def test_command_line_file_writes_keep_file_type_for_cromwell(self):
        expected = {"StageMethBatChromosomes": 1, "CompareMethBatBoundaries": 2}
        for task_name, number in expected.items():
            writes = []
            for part in self.tasks[task_name].command.parts:
                if not isinstance(part, WDL.Expr.Placeholder):
                    continue
                pending = [part.expr]
                while pending:
                    node = pending.pop()
                    if isinstance(node, WDL.Expr.Apply) and node.function_name == "write_lines":
                        # Terra maps the placeholder result to a local path only
                        # while the generated path still has File type.
                        self.assertIs(node, part.expr, task_name)
                        self.assertIsInstance(part.expr.type, WDL.Type.File, task_name)
                        writes.append(node)
                    pending.extend(node.children)
            self.assertEqual(len(writes), number, task_name)

    def test_disjoint_label_balanced_split_and_missing_localization(self):
        manifest = self.root / "manifest with 'quote' $ literal.tsv"
        rows = [f"id{i:03d}\tgs://test-data/sample{i:03d}\t{'F' if i < 115 else 'M'}" for i in range(231)]
        manifest.write_text("identifier\tfilename\tlabels\n" + "\n".join(rows) + "\n")
        uri = "gs://test-data/cohort.tsv"
        inputs = {"manifest_file": uri, "subset_count": 2, "subset_sizes": [], "split_seed": 42}
        first, second = self.root / "first", self.root / "second"
        for directory in (first, second):
            self.run_task("SplitMethBatCohort", directory, inputs, {uri: str(manifest)})
        manifests = sorted((first / "subsets").glob("*.tsv"))
        groups = [path.read_text().splitlines()[1:] for path in manifests]
        self.assertEqual(sorted(map(len, groups)), [115, 116])
        self.assertEqual(sorted(groups[0] + groups[1]), rows)
        self.assertTrue(set(groups[0]).isdisjoint(groups[1]))
        self.assertEqual(sorted(sum(row.endswith("\tF") for row in group) for group in groups), [57, 58])
        for path in manifests:
            self.assertEqual(path.read_bytes(), (second / "subsets" / path.name).read_bytes())
        self.run_task("SplitMethBatCohort", self.root / "bad", inputs, {}, "localization error")
        pilot = self.root / "pilot"
        pilot_inputs = dict(inputs, subset_sizes=[25, 50, 100])
        self.run_task("SplitMethBatCohort", pilot, pilot_inputs, {uri: str(manifest)})
        groups = [path.read_text().splitlines()[1:] for path in sorted((pilot / "subsets").glob("*.tsv"))]
        self.assertEqual([len(group) for group in groups], [100, 100])
        self.assertTrue(set(groups[0]).isdisjoint(groups[1]))
        for size in (25, 50, 100):
            for group in groups:
                n_f = sum(row.endswith("\tF") for row in group[:size])
                self.assertLessEqual(abs(n_f - (size - n_f)), 1)
        self.run_task(
            "SplitMethBatCohort", self.root / "too_large", dict(inputs, subset_sizes=[116]),
            {uri: str(manifest)}, "232 unique samples",
        )
        manifest.write_text("id000\t\tF\n")
        self.run_task("SplitMethBatCohort", self.root / "bad_row", inputs, {uri: str(manifest)}, "pb-CpG-tools")

    def install_stub(self, name, source):
        path = self.bin / name
        path.write_text("#!/usr/bin/env python3\n" + source)
        path.chmod(0o755)

    def test_stage_and_segment_chromosome_with_cloud_paths(self):
        self.install_stub("gsutil", """import pathlib, shutil, sys
args = sys.argv[1:]
assert args[:4] == ['-q', '-u', 'test-project', 'cp'], args
for uri in args[4:-1]:
    assert uri.startswith('gs://test-data/') and '*' not in uri, uri
    name = uri.rsplit('/', 1)[1]
    shutil.copyfile(pathlib.Path('../objects', name), pathlib.Path(args[-1], name))
""")
        if not os.environ.get("METHBAT_REAL"):
            self.install_stub("methbat", """import gzip, pathlib, sys
args = sys.argv[1:]
if args == ['--version']:
    print('methbat 1.2.0-7a582da')
    sys.exit(0)
assert args[:3] == ['joint-segment', '-t', '2']
assert '--enable-nodata-segments' in args and '--condense-bed-labels' not in args
assert args[args.index('--min-sample-frac') + 1] == '0.9'
assert args[args.index('--min-asm-frac') + 1] == '0.6'
manifest = pathlib.Path(args[args.index('--input-collection') + 1]).read_text().splitlines()
assert len(manifest) == 2
_, pileup, label = manifest[1].split('\\t')
assert label in ('F;study', '')
rows = gzip.open(pileup, 'rt').read().splitlines()
assert rows[0] == '##base_modification=5mC'
assert rows[1].split('\\t')[7] == 'type'
assert {row.split('\\t')[7] for row in rows[2:]} == ({'Total', 'hap1'} if label else {'Total'})
assert all(row.startswith('chr21\\t') for row in rows[2:])
out = args[args.index('--output-prefix') + 1]
for suffix in ('meth_regions.bed', 'combined_methyl.bedgraph', 'asm.bedgraph'):
    pathlib.Path(f'{out}.{suffix}').write_text('chr21\\t10\\t11\\t1\\n')
""")
        objects = self.root / "objects"
        objects.mkdir()
        for track in ("combined", "hap1", "hap2"):
            path = objects / f"sample.{track}.bed.gz"
            header = "#chrom\tstart\tend\tmod_score\ttype\tcov\tmod_count\tunmod_count\tavg_mod_score\tavg_unmod_score\n"
            payload = header + f"chr20\t20\t21\t10.0\t{track}\t10\t1\t9\t0.9\t0.1\n"
            if track != "hap2":
                payload += "".join(
                    f"chr21\t{pos}\t{pos+1}\t10.0\t{track}\t10\t1\t9\t0.9\t0.1\n"
                    for pos in range(10, 2010, 10)
                )
            with path.open("wb") as handle:
                subprocess.run(["bgzip", "-c"], input=payload.encode(), stdout=handle, check=True)
            subprocess.run(["tabix", "-p", "bed", path], check=True)
        manifest = self.root / "subset local.tsv"
        for path in sorted(objects.glob("sample.*")):
            (objects / path.name.replace("sample.", "other.", 1)).write_bytes(path.read_bytes())
        manifest.write_text(
            "identifier\tfilename\tlabels\nA\tgs://test-data/sample\tF;study\nB\tgs://test-data/other\tM\n"
        )
        staging = self.root / "stage"
        subset_uri = "gs://execution-outputs/subset-000001.tsv"
        stage_inputs = {
            "subset_manifest": subset_uri, "chromosomes": ["chr21", "chr20"],
            "transfer_workers": 2, "gcs_billing_project": "test-project",
        }
        self.run_task("StageMethBatChromosomes", staging, stage_inputs, {subset_uri: str(manifest)})
        bundle = staging / "bundles/chr21.tar"
        with tarfile.open(bundle) as archive:
            self.assertIn("collection.tsv", archive.getnames())
        self.assertTrue((staging / "bundles/chr20.tar").is_file())
        self.run_task("StageMethBatChromosomes", self.root / "bad_stage", stage_inputs, {}, "localization error")
        self.run_task("StageMethBatChromosomes", self.root / "absent", dict(stage_inputs, chromosomes=["chr99"]), {subset_uri: str(manifest)}, "No combined CpG rows")
        self.run_task("StageMethBatChromosomes", self.root / "invalid", dict(stage_inputs, chromosomes=["chr21\nchr20"]), {subset_uri: str(manifest)}, "Invalid chromosome")
        uri = "gs://execution-outputs/chr21.tar"
        segment_inputs = {
            "chromosome_bundle": uri, "subset_name": "subset-000001", "chromosome": "chr21",
            "threads": 2, "sample_limit": 1, "condense_bed_labels": False, "enable_nodata_segments": True,
            "docker_image": "test-methbat",
        }
        segmentation = self.root / "segment"
        self.run_task("SegmentMethBatChromosome", segmentation, segment_inputs, {uri: str(bundle)})
        prefix = segmentation / "output/subset-000001.chr21.n1"
        for suffix in ("meth_regions.bed", "aggregate.meth_regions.bed", "aggregate.5mC.bed.gz",
                       "aggregation_qc.tsv", "method_comparison.tsv", "method_differences.tsv"):
            self.assertTrue(Path(f"{prefix}.{suffix}").is_file(), suffix)
        with Path(f"{prefix}.method_comparison.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["exact_rows_match"] == "true" for row in rows), rows)
        self.run_task("SegmentMethBatChromosome", self.root / "bad_segment", segment_inputs, {}, "localization error")

    def test_aggregation_applies_both_cutoffs_and_uses_equal_donor_weights(self):
        work = self.root / "aggregation"
        tracks = work / "tracks"
        tracks.mkdir(parents=True)
        samples = []
        # Sites 10 and 20 fall on either side of the 90% combined cutoff.
        # Sites 20 and 30 fall on either side of the 60% paired cutoff.
        for donor in range(25):
            prefix = f"tracks/sample-{donor+1:06d}"
            samples.append(f"donor{donor}\t{prefix}\tX\n")
            for i, track in enumerate(("combined", "hap1", "hap2")):
                with gzip.open(work / f"{prefix}.{track}.bed.gz", "wt") as handle:
                    for start, n_combined, n_paired in ((10, 22, 15), (20, 23, 15), (30, 25, 14)):
                        if i == 0 and donor < n_combined:
                            # The higher-depth call must still have equal weight.
                            counts = (0, 100) if donor == 0 else (1, 2)
                        elif i in (1, 2) and donor < n_paired:
                            counts = (1, 1) if (donor + i) % 2 else (0, 1)
                        elif i == 1 and donor == n_paired:
                            # An unpaired call adds no haplotype support.
                            counts = (1, 1)
                        else:
                            continue
                        modified, depth = counts
                        handle.write(f"chr21\t{start}\t{start+1}\t0.0\t{track}\t{depth}\t{modified}"
                                     f"\t{depth-modified}\t0.9\t0.1\n")
        (work / "collection.tsv").write_text("identifier\tfilename\tlabels\n" + "".join(samples))
        subprocess.run(
            ["methbat_prepare_pileups.py", "--collection", "collection.tsv", "--chromosome", "chr21",
             "--sample-limit", "25", "--output-prefix", "pilot"], cwd=work, env=self.env, check=True,
        )
        with gzip.open(work / "pilot.aggregate.5mC.bed.gz", "rt") as handle:
            rows = {(fields[1], fields[7]): fields[8:11]
                    for fields in (line.split("\t") for line in handle) if fields[0] == "chr21"}
        self.assertEqual(rows, {
            ("20", "Total"): ["10000", "4783", "5217"],
            ("20", "hap1"): ["10000", "0", "10000"],
            ("20", "hap2"): ["10000", "10000", "0"],
            ("30", "Total"): ["10000", "4800", "5200"],
        })
        with (work / "pilot.aggregation_qc.tsv").open() as handle:
            qc = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([qc[field] for field in ("sample_count", "candidate_sites", "combined_pass_sites",
                                                  "paired_haplotype_pass_sites")], ["25", "3", "2", "1"])

    @unittest.skipUnless(os.environ.get("METHBAT_REAL"), "requires the pinned MethBat binary")
    def test_real_methbat_on_25_donors(self):
        fixture = self.root / "donors"
        (fixture / "tracks").mkdir(parents=True)
        samples = []
        for donor in range(25):
            prefix = f"tracks/sample-{donor+1:06d}"
            samples.append(f"donor{donor}\t{prefix}\tX\n")
            for index, track in enumerate(("combined", "hap1", "hap2")):
                with gzip.open(fixture / f"{prefix}.{track}.bed.gz", "wt") as handle:
                    blocks = ((25, 25, 90), (22, 15, 10), (23, 15, 40), (25, 14, 10))
                    for block, (n_total, n_paired, combined_pct) in enumerate(blocks, 1):
                        if donor >= (n_total if index == 0 else n_paired):
                            continue
                        depth = 100 * (donor % 3 + 1)
                        percent = combined_pct if index == 0 else (90 if (index + donor) % 2 else 10)
                        modified = percent * depth // 100
                        for start in range(block * 3000, block * 3000 + 1500, 10):
                            handle.write(f"chr21\t{start}\t{start+1}\t{percent}.0\t{track}\t{depth}"
                                         f"\t{modified}\t{depth-modified}\t0.9\t0.1\n")
        (fixture / "collection.tsv").write_text("identifier\tfilename\tlabels\n" + "".join(samples))
        bundle = self.root / "chr21.tar"
        with tarfile.open(bundle, "w") as archive:
            archive.add(fixture / "collection.tsv", "collection.tsv")
            archive.add(fixture / "tracks", "tracks")
        uri = "gs://execution-outputs/chr21.tar"
        work = self.root / "native"
        self.run_task("SegmentMethBatChromosome", work, {
            "chromosome_bundle": uri, "subset_name": "subset-000001", "chromosome": "chr21",
            "threads": 2, "sample_limit": 25, "condense_bed_labels": False, "enable_nodata_segments": True,
            "docker_image": "test-methbat",
        }, {uri: str(bundle)})
        prefix = work / "output/subset-000001.chr21.n25"
        with Path(f"{prefix}.method_comparison.tsv").open() as handle:
            summary = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(row["exact_rows_match"] == "true" for row in summary), summary)
        self.assertTrue(all(int(row["joint_rows"]) > 0 for row in summary), summary)

    def test_method_report_detects_label_boundary_and_score_changes(self):
        work = self.root / "differences"
        work.mkdir()
        files = {
            "meth_regions.bed": ("chr21\t10\t30\tM\n", "chr21\t10\t30\tU\n"),
            "combined_methyl.bedgraph": ("chr21\t10\t30\t51.0\n", "chr21\t11\t30\t51.0\n"),
            "asm.bedgraph": ("chr21\t10\t30\t50.0\n", "chr21\t10\t30\t49.9\n"),
        }
        for suffix, (joint, aggregate) in files.items():
            (work / f"pilot.{suffix}").write_text(joint)
            (work / f"pilot.aggregate.{suffix}").write_text(aggregate)
        subprocess.run(["methbat_compare_methods.py", "--prefix", "pilot", "--subset", "subset-000001",
                        "--chromosome", "chr21", "--sample-count", "25"], cwd=work, env=self.env, check=True)
        with (work / "pilot.method_comparison.tsv").open() as handle:
            summary = {row["output"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertTrue(all(row["exact_rows_match"] == "false" for row in summary.values()))
        self.assertEqual(summary["meth_regions.bed"]["joint_only_rows"], "1")
        self.assertEqual(summary["combined_methyl.bedgraph"]["shared_coordinates"], "0")
        self.assertEqual(summary["asm.bedgraph"]["max_abs_score_delta_pp"], "0.1")
        self.assertEqual(len((work / "pilot.method_differences.tsv").read_text().splitlines()), 7)

    def test_compare_boundaries_by_size_and_label(self):
        beds = {
            "25a": "chr21\t10\t30\tM\nchr21\t40\t60\tM\nchr21\t100\t120\tASM\n",
            "25b": "chr21\t12\t27\tM\nchr21\t55\t63\tM\nchr21\t100\t120\tU\n",
            "50a": "chr21\t10\t30\tM\n",
            "50b": "chr21\t10\t30\tM\n",
            "100a": "",
            "100b": "# No calls\n",
        }
        keys = ["subset-000001\tchr21\t25\t25", "subset-000002\tchr21\t25\t25",
                "subset-000001\tchr21\t50\t50", "subset-000002\tchr21\t50\t50",
                "subset-000001\tchr21\t100\t100", "subset-000002\tchr21\t100\t100"]
        uris = {f"gs://execution-outputs/{name}.bed": self.root / f"{name} 'local'.bed" for name in beds}
        for path, content in zip(uris.values(), beds.values()):
            path.write_text(content)
        inputs = {"region_beds": list(uris), "run_keys": keys, "min_reciprocal_overlap": 0.5}
        output = self.root / "comparison"
        self.run_task("CompareMethBatBoundaries", output, inputs, {uri: str(path) for uri, path in uris.items()})
        with (output / "boundary_summary.tsv").open() as handle:
            summary = {(row["subset_size"], row["label"]): row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(summary["25", "M"]["matched_regions"], "1")
        self.assertEqual(summary["25", "M"]["unmatched_a"], "1")
        self.assertEqual(summary["25", "M"]["unmatched_b"], "1")
        self.assertEqual(summary["25", "M"]["median_abs_end_shift_bp"], "3")
        self.assertEqual(summary["25", "ASM"]["matched_regions"], "0")
        self.assertEqual(summary["50", "M"]["median_abs_start_shift_bp"], "0")
        self.assertEqual(summary["100", ""]["regions_a"], "0")
        self.assertEqual(summary["100", ""]["median_abs_start_shift_bp"], "NA")
        with (output / "boundary_shifts.tsv").open() as handle:
            shifts = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([(row["start_shift_bp"], row["end_shift_bp"]) for row in shifts], [("2", "-3"), ("0", "0")])
        self.run_task("CompareMethBatBoundaries", self.root / "bad_compare", inputs, {}, "localization error")


if __name__ == "__main__":
    unittest.main()
