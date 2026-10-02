"""Check HiRes gene chunks and merged estimates with synthetic task inputs.

File localization and command-time File writers are modeled explicitly.
The native calculation is a process fixture; no credentials or Terra job are used.
"""
import csv
import hashlib
import io
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]
WDL_PATH = ROOT / "workflows/cell_type_specific_expression/cibersortx_hires.wdl"
GENES = ["G3", "G1", "G4", "G2", "G5"]
CELLS = ["B", "CD4_T"]
SAMPLES = ["s%d" % n for n in range(1, 9)]
VALUES = ["1", "NA", "NaN", "", "1.000", "2e-03", ".400", "-0.5"]


class CloudGeneratedStdLib(WDL.StdLib.Base):
    def __init__(self, directory):
        super().__init__("1.0", write_dir=str(directory))
        self.generated_paths = {}

    def _virtualize_filename(self, filename):
        uri = "gs://test-bucket/generated/" + Path(filename).name
        self.generated_paths[uri] = filename
        return uri

    def _devirtualize_filename(self, filename):
        return self.generated_paths.get(filename, filename)


class ChunkingTest(unittest.TestCase):
    def task(self, name):
        document = WDL.load(str(WDL_PATH))
        tasks = {task.name: task for task in document.tasks}
        self.assertIn(name, tasks, "The chunking task has not been created: " + name)
        return tasks[name]

    def execute(self, name, base, inputs, native=None):
        """Render after localization and map command-time generated File values.

        Files inside a previously serialized list are never rewritten. Thus a
        task declaration that writes cloud paths fails this same boundary test.
        """
        base.mkdir(parents=True, exist_ok=True)
        task = self.task(name)
        env = WDL.values_from_json(inputs, task.available_inputs)
        downloads = {}
        def cloud(value):
            if "://" in value.value:
                return value.value
            uri = "gs://test-bucket/inputs/%d/%s" % (len(downloads), Path(value.value).name)
            downloads[uri] = value.value
            return uri
        env = WDL.Value.rewrite_env_paths(env, cloud)
        stdlib = CloudGeneratedStdLib(base / "generated")
        for declaration in task.inputs:
            if not env.has_binding(declaration.name):
                if declaration.expr is not None:
                    env = env.bind(declaration.name, declaration.expr.eval(env, stdlib))
                elif declaration.type.optional:
                    env = env.bind(declaration.name, WDL.Value.Null())
        for declaration in task.postinputs:
            env = env.bind(declaration.name, declaration.expr.eval(env, stdlib))
        env = WDL.Value.rewrite_env_paths(env, lambda value: downloads.get(
            value.value, stdlib.generated_paths.get(value.value, value.value)))
        rendered = []
        # Cromwell maps a File returned by a command placeholder. It cannot map
        # a cloud path after the placeholder has converted it to a String.
        for part in WDL.Expr.String._dedent(task.command.parts):
            if not isinstance(part, WDL.Expr.Placeholder):
                rendered.append(part)
                continue
            value = part.expr.eval(env, stdlib)
            if isinstance(value, WDL.Value.File):
                self.assertFalse(part.options, "A generated File placeholder must use its local path")
                localized = WDL.Value.rewrite_paths(value, lambda file: downloads.get(
                    file.value, stdlib.generated_paths.get(file.value, file.value)))
                rendered.append(localized.value)
            else:
                if "write_lines" in str(part.expr):
                    self.fail("A command-time file list was converted to String before localization")
                rendered.append(part.eval(env, stdlib).value)
        command = "".join(rendered)
        src = base / "src"
        src.mkdir(exist_ok=True)
        executable = src / "CIBERSORTxHiRes"
        if native is not None:
            executable.write_text(native)
            executable.chmod(0o755)
        image = os.environ.get("CIBERSORTX_TEST_IMAGE")
        if not image:
            command = command.replace("/src", str(src))
        script = base / "command.sh"
        script.write_text(command)
        invocation = ["bash", str(script)]
        if image:
            invocation = ["docker", "run", "--rm", "--platform", "linux/amd64",
                          "--entrypoint", "/bin/bash", "--workdir", str(base),
                          "-v", str(base.parent) + ":" + str(base.parent)]
            if native is not None:
                invocation.extend(["-v", str(executable) + ":/src/CIBERSORTxHiRes:ro"])
            invocation.extend([image, "-c",
                               'set +e; bash "$1"; result=$?; chown -R "$2:$3" "$4"; exit "$result"',
                               "chunk-test", str(script), str(os.getuid()), str(os.getgid()), str(base.parent)])
        process = subprocess.run(invocation, cwd=base, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
        return process, stdlib

    @staticmethod
    def panel_file(base, text=None):
        base.mkdir(parents=True, exist_ok=True)
        file = base / "panel's $(touch INJECTED).txt"
        file.write_text(text if text is not None else "\n".join(GENES) + "\n")
        return file

    def split(self, base, text=None, size=2):
        panel = self.panel_file(base / "source", text)
        process, _ = self.execute("SplitGeneSubset", base / "split", {
            "gene_subset": str(panel), "genes_per_chunk": size})
        chunks = sorted((base / "split/chunks").glob("genes_*.txt"))
        return process, panel, chunks

    def test_split_contiguous_chunks_keeps_order_and_last_short_chunk(self):
        with tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
            base = Path(temp)
            process, _, chunks = self.split(base)
            self.assertEqual(process.returncode, 0, process.stdout)
            self.assertEqual([chunk.read_text().splitlines() for chunk in chunks],
                             [["G3", "G1"], ["G4", "G2"], ["G5"]])
            self.assertTrue(all(0 < len(chunk.read_text().splitlines()) <= 2 for chunk in chunks))
            self.assertTrue((base / "split/chunk_plan.json").is_file())
            json.loads((base / "split/chunk_plan.json").read_text())
            self.assertTrue((base / "split/split.log").is_file())
            self.assertFalse(list(base.rglob("INJECTED")))

    def test_split_small_panel_produces_one_nonempty_chunk(self):
        with tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
            process, _, chunks = self.split(Path(temp), size=1000)
            self.assertEqual(process.returncode, 0, process.stdout)
            self.assertEqual([chunk.read_text().splitlines() for chunk in chunks], [GENES])

    def test_split_size_one_and_exactly_divisible_panel_have_no_empty_chunk(self):
        cases = [(GENES, 1, [[gene] for gene in GENES]),
                 (GENES[:4], 2, [GENES[:2], GENES[2:4]])]
        for genes, size, expected in cases:
            with self.subTest(size=size), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                process, _, chunks = self.split(Path(temp), "\n".join(genes) + "\n", size)
                self.assertEqual(process.returncode, 0, process.stdout)
                self.assertEqual([chunk.read_text().splitlines() for chunk in chunks], expected)

    def test_split_invalid_panel_or_size_fails(self):
        cases = [("", 2), ("G1\nG1\n", 2), ("G1\n\nG2\n", 2),
                 (" G1\n", 2), ("G1\tG2\n", 2), ("G1\n", 0), ("G1\n", -1)]
        for text, size in cases:
            with self.subTest(text=text, size=size), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                process, _, _ = self.split(Path(temp), text=text, size=size)
                self.assertNotEqual(process.returncode, 0, process.stdout)

    @staticmethod
    def matrix_text(genes, values=None):
        values = VALUES if values is None else values
        return "GeneSymbol\t" + "\t".join(SAMPLES) + "\n" + "".join(
            gene + "\t" + "\t".join(values) + "\n" for gene in genes)

    def merge_fixture(self, base):
        panel = self.panel_file(base / "source")
        inputs = {"gene_subset": str(panel), "chunk_gene_subsets": [],
                  "chunk_expression_matrices": [], "input_validations": [], "run_logs": [],
                  "memory_gb": 4, "disk_gb": 20}
        for index, genes in enumerate((["G3", "G1"], ["G4", "G2"], ["G5"])):
            directory = base / "localized chunk's" / ("chunk%d" % index)
            directory.mkdir(parents=True)
            chunk = directory / "genes.txt"
            chunk.write_text("\n".join(genes) + "\n")
            inputs["chunk_gene_subsets"].append(str(chunk))
            report = {
                "sample_count": len(SAMPLES), "samples": SAMPLES, "cell_types": CELLS,
                "mixture_gene_count": 6, "signature_gene_count": 1, "subset_genes": genes,
                "threads": 8, "nsampling": 1, "nsampling2": 1,
                "fraction_input_sample_count": 8, "fraction_retained_sample_count": 8,
                "fraction_dropped_sample_count": 0,
                "negative_mixture_values": 0, "negative_signature_values": 0,
                "sha256": {"mixture": "0" * 64, "signature": "1" * 64, "fractions": "2" * 64,
                           "gene_subset": hashlib.sha256(chunk.read_bytes()).hexdigest()},
            }
            validation = directory / "input_validation.json"
            validation.write_text(json.dumps(report))
            inputs["input_validations"].append(str(validation))
            run_log = directory / "hires.log"
            run_log.write_text("Synthetic chunk %d completed.\n" % index)
            inputs["run_logs"].append(str(run_log))
            for cell in CELLS:
                matrix = directory / ("CIBERSORTxHiRes_ZNF804A_controls_%s_Window8.txt" % cell)
                matrix.write_text(self.matrix_text(list(reversed(genes))))
                inputs["chunk_expression_matrices"].append(str(matrix))
        return inputs

    def merged_files(self, base):
        return sorted((base / "results").glob("CIBERSORTxHiRes*_Window*.txt"))

    def test_merge_restores_panel_order_and_preserves_numeric_and_missing_text(self):
        with tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
            base = Path(temp)
            inputs = self.merge_fixture(base)
            # Equal basenames must not overwrite matrices from another chunk.
            self.assertEqual(len({Path(path).name for path in inputs["chunk_expression_matrices"]}), 2)
            process, stdlib = self.execute("MergeHiRes", base / "merge", inputs)
            self.assertEqual(process.returncode, 0, process.stdout)
            files = self.merged_files(base / "merge")
            self.assertEqual(len(files), len(CELLS))
            for file in files:
                self.assertEqual(file.read_text(), self.matrix_text(GENES))
            self.assertTrue(stdlib.generated_paths, "Merge array File inputs must exercise command-time lists")
            for path in stdlib.generated_paths.values():
                self.assertNotIn("gs://", Path(path).read_text(), "Array File values were serialized before localization")
            report = json.loads((base / "merge/input_validation.json").read_text())
            self.assertEqual(report["subset_genes"], GENES)
            self.assertEqual(report["samples"], SAMPLES)
            self.assertEqual(report["cell_types"], CELLS)
            self.assertTrue((base / "merge/output_validation.json").is_file())
            self.assertTrue((base / "merge/hires.log").is_file())

    def test_merge_rejects_wrong_chunk_partitions(self):
        for problem in ("overlap", "missing", "extra", "duplicate_panel"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                base = Path(temp)
                inputs = self.merge_fixture(base)
                if problem == "overlap":
                    Path(inputs["chunk_gene_subsets"][1]).write_text("G1\nG4\nG2\n")
                elif problem == "missing":
                    for name in ("chunk_gene_subsets", "input_validations", "run_logs"):
                        inputs[name].pop()
                    inputs["chunk_expression_matrices"] = inputs["chunk_expression_matrices"][:-2]
                elif problem == "extra":
                    Path(inputs["chunk_gene_subsets"][-1]).write_text("G5\nEXTRA\n")
                else:
                    Path(inputs["gene_subset"]).write_text("\n".join(GENES + ["G3"]) + "\n")
                process, _ = self.execute("MergeHiRes", base / "merge", inputs)
                self.assertNotEqual(process.returncode, 0, process.stdout)

    def test_merge_rejects_changed_samples_cells_windows_or_missing_matrices(self):
        for problem in ("samples", "cell", "window", "missing_matrix", "duplicate_matrix"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                base = Path(temp)
                inputs = self.merge_fixture(base)
                file = Path(inputs["chunk_expression_matrices"][2])
                if problem == "samples":
                    file.write_text(file.read_text().replace("s8", "WRONG", 1))
                elif problem in ("cell", "window"):
                    renamed = file.with_name(file.name.replace("_B_", "_NK_") if problem == "cell"
                                             else file.name.replace("Window8", "Window9"))
                    file.rename(renamed)
                    inputs["chunk_expression_matrices"][2] = str(renamed)
                elif problem == "missing_matrix":
                    inputs["chunk_expression_matrices"].pop()
                else:
                    inputs["chunk_expression_matrices"][3] = inputs["chunk_expression_matrices"][2]
                process, _ = self.execute("MergeHiRes", base / "merge", inputs)
                self.assertNotEqual(process.returncode, 0, process.stdout)

    def test_merge_rejects_duplicate_missing_extra_or_nonrectangular_matrix_rows(self):
        for problem in ("duplicate", "missing", "extra", "nonrectangular"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                base = Path(temp)
                inputs = self.merge_fixture(base)
                path = Path(inputs["chunk_expression_matrices"][0])
                lines = path.read_text().splitlines()
                if problem == "duplicate":
                    lines.append(lines[1])
                elif problem == "missing":
                    lines.pop()
                elif problem == "extra":
                    lines.append("EXTRA\t" + "\t".join(VALUES))
                else:
                    lines[1] = lines[1].rsplit("\t", 1)[0]
                path.write_text("\n".join(lines) + "\n")
                process, _ = self.execute("MergeHiRes", base / "merge", inputs)
                self.assertNotEqual(process.returncode, 0, process.stdout)

    def test_merge_rejects_nonnumeric_or_nonfinite_estimates(self):
        for value in ("not_numeric", "inf", "-Infinity"):
            with self.subTest(value=value), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                base = Path(temp)
                inputs = self.merge_fixture(base)
                path = Path(inputs["chunk_expression_matrices"][0])
                lines = path.read_text().splitlines()
                row = lines[1].split("\t")
                row[1] = value
                lines[1] = "\t".join(row)
                path.write_text("\n".join(lines) + "\n")
                process, _ = self.execute("MergeHiRes", base / "merge", inputs)
                self.assertNotEqual(process.returncode, 0, process.stdout)
                self.assertRegex(process.stdout, r"(?i)(nonnumeric|nonfinite)")

    def test_reversed_chunk_arrays_still_restore_original_panel_order(self):
        with tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
            base = Path(temp)
            inputs = self.merge_fixture(base)
            for name in ("chunk_gene_subsets", "input_validations", "run_logs"):
                inputs[name].reverse()
            matrices = inputs["chunk_expression_matrices"]
            inputs["chunk_expression_matrices"] = matrices[4:6] + matrices[2:4] + matrices[0:2]
            process, _ = self.execute("MergeHiRes", base / "merge", inputs)
            self.assertEqual(process.returncode, 0, process.stdout)
            for file in self.merged_files(base / "merge"):
                self.assertEqual(file.read_text(), self.matrix_text(GENES))

    def test_merge_rejects_changed_background_hashes_settings_or_subset_reports(self):
        for problem in ("mixture", "signature", "fractions", "nsampling", "subset_genes", "subset_hash"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                base = Path(temp)
                inputs = self.merge_fixture(base)
                path = Path(inputs["input_validations"][1])
                report = json.loads(path.read_text())
                if problem in ("mixture", "signature", "fractions"):
                    report["sha256"][problem] = "f" * 64
                elif problem == "subset_hash":
                    report["sha256"]["gene_subset"] = "f" * 64
                elif problem == "nsampling":
                    report[problem] = 100
                else:
                    report[problem] = ["G4", "WRONG"]
                path.write_text(json.dumps(report))
                process, _ = self.execute("MergeHiRes", base / "merge", inputs)
                self.assertNotEqual(process.returncode, 0, process.stdout)

    def test_unresolved_cloud_inputs_fail_split_and_merge(self):
        with tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
            base = Path(temp)
            process, _ = self.execute("SplitGeneSubset", base / "split", {
                "gene_subset": "gs://bucket/panel.txt", "genes_per_chunk": 2})
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("Localization error", process.stdout)
        for name in ("gene_subset", "chunk_gene_subsets", "chunk_expression_matrices", "input_validations", "run_logs"):
            with self.subTest(file=name), tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
                base = Path(temp)
                inputs = self.merge_fixture(base)
                if name == "gene_subset":
                    inputs[name] = "gs://bucket/unresolved.txt"
                else:
                    inputs[name][0] = "gs://bucket/unresolved.txt"
                process, _ = self.execute("MergeHiRes", base / "merge", inputs)
                self.assertNotEqual(process.returncode, 0)
                self.assertIn("Localization error", process.stdout)

    @staticmethod
    def native_fixture():
        return "#!/usr/bin/env python3\n" + "\n".join([
            "import hashlib, json, pathlib, sys",
            "args = dict(zip(sys.argv[1::2], sys.argv[2::2]))",
            "out = pathlib.Path(args['--outdir']).resolve()",
            "data = pathlib.Path.cwd() / 'data'",
            "backgrounds = {}",
            "for name, flag in [('mixture', '--mixture'), ('signature', '--sigmatrix')]:",
            "    backgrounds[name] = hashlib.sha256((data / args[flag]).read_bytes()).hexdigest()",
            "backgrounds['fractions'] = hashlib.sha256((out / args['--cibresults']).read_bytes()).hexdigest()",
            "(out.parent / 'native_backgrounds.json').write_text(json.dumps(backgrounds))",
            "genes = (data / args['--subsetgenes']).read_text().splitlines()",
            "header = 'GeneSymbol\\t' + '\\t'.join(" + repr(SAMPLES) + ") + '\\n'",
            "rows = ''.join(gene + '\\t' + '\\t'.join(" + repr(VALUES) + ") + '\\n' for gene in reversed(genes))",
            "for cell in " + repr(CELLS) + ":",
            "    (out / ('CIBERSORTxHiRes_' + args['--label'] + '_' + cell + '_Window8.txt')).write_text(header + rows)",
            "print('Synthetic native chunk completed.', flush=True)",
            "",
        ])

    def test_synthetic_task_chain_keeps_full_backgrounds_for_every_chunk(self):
        with tempfile.TemporaryDirectory(prefix="hires-chunks-") as temp:
            base = Path(temp)
            process, panel, chunks = self.split(base)
            self.assertEqual(process.returncode, 0, process.stdout)
            source = base / "source"
            mixture = source / "mixture.txt"
            mixture.write_text("GeneSymbol\t" + "\t".join(SAMPLES) + "\n" +
                               "".join(gene + "\t" + "\t".join(["2"] * 8) + "\n" for gene in ["sig"] + GENES))
            signature = source / "signature.txt"
            signature.write_text("GeneSymbol\tB\tCD4_T\nsig\t1\t2\n")
            fractions = source / "fractions.txt"
            fractions.write_text("Mixture\tB\tCD4_T\n" + "".join(sample + "\t0.25\t0.75\n" for sample in SAMPLES))
            shared = {"mixture": str(mixture), "signature": str(signature), "fractions": str(fractions),
                      "username": "fixture@example.org", "token": "fixture-token", "threads": 8,
                      "nsampling": 1, "nsampling2": 1, "memory_gb": 16, "disk_gb": 20}
            expected_hashes = {name: hashlib.sha256(Path(shared[name]).read_bytes()).hexdigest()
                               for name in ("mixture", "signature", "fractions")}
            merge = {"gene_subset": str(panel), "chunk_gene_subsets": [str(chunk) for chunk in chunks],
                     "chunk_expression_matrices": [], "input_validations": [], "run_logs": [], "memory_gb": 4, "disk_gb": 20}
            for index, chunk in enumerate(chunks):
                task_root = base / ("native_chunk%d" % index)
                process, _ = self.execute("RunHiRes", task_root, dict(shared, gene_subset=str(chunk)), self.native_fixture())
                self.assertEqual(process.returncode, 0, process.stdout)
                self.assertEqual(json.loads((task_root / "native_backgrounds.json").read_text()), expected_hashes)
                report = json.loads((task_root / "input_validation.json").read_text())
                self.assertEqual({name: report["sha256"][name] for name in expected_hashes}, expected_hashes)
                matrices = sorted((task_root / "results").glob("CIBERSORTxHiRes*_Window*.txt"))
                # Cell order, rather than sorted file names, defines each block.
                merge["chunk_expression_matrices"].extend(str(next(file for file in matrices
                    if "_" + cell + "_Window" in file.name)) for cell in CELLS)
                merge["input_validations"].append(str(task_root / "input_validation.json"))
                merge["run_logs"].append(str(task_root / "hires.log"))
            process, _ = self.execute("MergeHiRes", base / "merge", merge)
            self.assertEqual(process.returncode, 0, process.stdout)
            files = self.merged_files(base / "merge")
            self.assertEqual(len(files), len(CELLS))
            for file in files:
                self.assertEqual(file.read_text(), self.matrix_text(GENES))
            for name, digest in expected_hashes.items():
                self.assertEqual(hashlib.sha256(Path(shared[name]).read_bytes()).hexdigest(), digest)

    def test_workflow_call_expressions_forward_full_backgrounds_and_chunk_only_panel(self):
        # Evaluate the real call expressions. A correct task chain does not
        # prove that workflow wiring forwards the same inputs and settings.
        document = WDL.load(str(WDL_PATH))
        workflow = document.workflow
        values = {
            "mixture": "gs://test-bucket/full_mixture.txt",
            "signature": "gs://test-bucket/full_signature.txt",
            "fractions": "gs://test-bucket/full_fractions.txt",
            "gene_subset": "gs://test-bucket/panel.txt",
            "username": "fixture@example.org", "token": "fixture-token",
            "threads": 6, "nsampling": 101, "nsampling2": 11,
            "memory_gb": 24, "disk_gb": 70, "genes_per_chunk": 2,
            "merge_memory_gb": 9, "merge_disk_gb": 55,
        }
        env = WDL.values_from_json(values, workflow.available_inputs)
        stdlib = WDL.StdLib.Base("1.0")
        split_call = next(node for node in workflow.body
                          if isinstance(node, WDL.Tree.Call) and node.callee.name == "SplitGeneSubset")
        split = {name: expression.eval(env, stdlib) for name, expression in split_call.inputs.items()}
        self.assertIsInstance(split["gene_subset"], WDL.Value.File)
        self.assertEqual(split["gene_subset"].value, values["gene_subset"])
        self.assertEqual(split["genes_per_chunk"].value, values["genes_per_chunk"])
        chunk_paths = ["gs://test-bucket/chunk%d/genes.txt" % index for index in range(3)]
        chunks = WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(path) for path in chunk_paths])
        env = env.bind(split_call.name + ".chunk_gene_subsets", chunks)
        scatter = next(node for node in workflow.body if isinstance(node, WDL.Tree.Scatter)
                       and any(isinstance(call, WDL.Tree.Call) and call.callee.name == "RunHiRes"
                               for call in node.body))
        run_call = next(call for call in scatter.body if isinstance(call, WDL.Tree.Call)
                        and call.callee.name == "RunHiRes")
        scattered = scatter.expr.eval(env, stdlib)
        self.assertEqual(scattered.json, chunk_paths)
        self.assertIsInstance(scattered.type.item_type, WDL.Type.File)
        for chunk in scattered.value:
            chunk_env = env.bind(scatter.variable, chunk)
            run = {name: expression.eval(chunk_env, stdlib) for name, expression in run_call.inputs.items()}
            self.assertEqual(run["gene_subset"].value, chunk.value)
            self.assertIsInstance(run["gene_subset"], WDL.Value.File)
            for name in ("mixture", "signature", "fractions"):
                self.assertIsInstance(run[name], WDL.Value.File)
                self.assertEqual(run[name].value, values[name])
            for name in ("username", "token", "threads", "nsampling", "nsampling2", "memory_gb", "disk_gb"):
                self.assertEqual(run[name].value, values[name])
        matrix_blocks = [["gs://test-bucket/chunk%d/%s.txt" % (index, cell) for cell in CELLS]
                         for index in range(3)]
        nested_matrices = WDL.Value.Array(WDL.Type.Array(WDL.Type.File()), [
            WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(path) for path in block])
            for block in matrix_blocks])
        validations = ["gs://test-bucket/chunk%d/input_validation.json" % index for index in range(3)]
        logs = ["gs://test-bucket/chunk%d/hires.log" % index for index in range(3)]
        env = env.bind(run_call.name + ".expression_matrices", nested_matrices)
        env = env.bind(run_call.name + ".input_validation", WDL.Value.Array(
            WDL.Type.File(), [WDL.Value.File(path) for path in validations]))
        env = env.bind(run_call.name + ".run_log", WDL.Value.Array(
            WDL.Type.File(), [WDL.Value.File(path) for path in logs]))
        merge_call = next(node for node in workflow.body
                          if isinstance(node, WDL.Tree.Call) and node.callee.name == "MergeHiRes")
        merge = {name: expression.eval(env, stdlib) for name, expression in merge_call.inputs.items()}
        self.assertIsInstance(merge["gene_subset"], WDL.Value.File)
        self.assertEqual(merge["gene_subset"].value, values["gene_subset"])
        self.assertEqual(merge["chunk_gene_subsets"].json, chunk_paths)
        self.assertEqual(merge["chunk_expression_matrices"].json, [path for block in matrix_blocks for path in block])
        self.assertEqual(merge["input_validations"].json, validations)
        self.assertEqual(merge["run_logs"].json, logs)
        for name in ("chunk_gene_subsets", "chunk_expression_matrices", "input_validations", "run_logs"):
            self.assertIsInstance(merge[name].type.item_type, WDL.Type.File)
        self.assertEqual(merge["memory_gb"].value, values["merge_memory_gb"])
        self.assertEqual(merge["disk_gb"].value, values["merge_disk_gb"])


if __name__ == "__main__":
    unittest.main()
