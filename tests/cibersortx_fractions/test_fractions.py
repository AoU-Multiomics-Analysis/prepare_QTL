"""Check rendered Fractions commands at file and process boundaries.

Cloud downloads are simulated by replacing File values with readable local
paths. The authenticated calculation is replaced with a native-process fixture.
These checks are not a Terra integration test or a scientific accuracy test.
"""
import csv
import importlib.util
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]
WDL_PATH = ROOT / "workflows/cell_type_specific_expression/cibersortx_fractions.wdl"
SAMPLES = ["s1", "s2", "s3"]
CELLS = ["B", "CD4_T"]


class FractionsTest(unittest.TestCase):
    def document(self):
        path = WDL_PATH
        self.assertTrue(path.is_file(), "The Terra Fractions WDL has not been created")
        return WDL.load(str(path))

    def run_task(self, change=None, native_status=0, output_problem=None, smode=False):
        doc = self.document()
        task = next(task for task in doc.tasks if task.name == "RunFractions")
        with tempfile.TemporaryDirectory(prefix="fractions localized ") as tmp:
            base = Path(tmp)
            localized = base / "localized input's"
            localized.mkdir()
            contents = {
                "mixture": "GeneSymbol\ts1\ts2\ts3\nsig\t1\t2\t3\nZNF804A\t2\t3\t4\n",
                "signature": "GeneSymbol\tB\tCD4_T\nsig\t1\t2\n",
                "refsample": "GeneSymbol\tB\tB.1\tCD4_T\tCD4_T.1\nsig\t1\t2\t3\t4\nZNF804A\t2\t3\t4\t5\n",
                "source_geps": "GeneSymbol\tB\tCD4_T\nsig\t1\t2\nZNF804A\t2\t3\n",
            }
            inputs = {}
            for name, content in contents.items():
                file = localized / (name + ".txt")
                file.write_text(content)
                if name in ("mixture", "signature") or smode:
                    inputs[name] = str(file)
            inputs.update(username="user'$(touch INJECTED)@example.org",
                          token="token'`touch INJECTED`", smode=smode,
                          quantile_normalization=False, permutations=100,
                          threads=8, memory_gb=16, disk_gb=30)
            if change:
                change(inputs, localized)
            # Start with cloud File values. Rewrite only those values to model
            # localization. Strings and absent optional Files stay unchanged.
            cloud_inputs = dict(inputs)
            downloads = {}
            for name in contents:
                declared = doc.workflow.available_inputs.resolve(name).type
                self.assertIsInstance(declared, WDL.Type.File)
                if name in inputs and inputs[name] is not None and "://" not in inputs[name]:
                    cloud_uri = "gs://test-bucket/" + name + ".txt"
                    downloads[cloud_uri] = inputs[name]
                    cloud_inputs[name] = cloud_uri
            env = WDL.values_from_json(cloud_inputs, task.available_inputs)
            env = WDL.Value.rewrite_env_paths(env, lambda file: downloads.get(file.value, file.value))
            stdlib = WDL.StdLib.Base("1.0")
            for decl in task.inputs:
                if not env.has_binding(decl.name):
                    if decl.expr is not None:
                        env = env.bind(decl.name, decl.expr.eval(env, stdlib))
                    elif decl.type.optional:
                        env = env.bind(decl.name, WDL.Value.Null())
            runtime_cpu = task.runtime["cpu"].eval(env, stdlib).value
            command = task.command.eval(env, stdlib).value
            src = base / "src"
            src.mkdir()
            native = src / "CIBERSORTxFractions"
            native.write_text(self.native_fixture(native_status, output_problem))
            native.chmod(0o755)
            # Map only the image's fixed /src path into the temporary test
            # filesystem. Keep the rendered command and File values intact.
            script = base / "command.sh"
            script.write_text(command.replace("/src", str(src)))
            proc = subprocess.run(["bash", str(script)], cwd=base, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
            args_file = base / "native_args.json"
            native_args = json.loads(args_file.read_text()) if args_file.exists() else None
            files = {str(file.relative_to(base)): file.read_text()
                     for file in base.rglob("*") if file.is_file()
                     and (file.suffix in (".txt", ".json", ".log"))
                     and "localized input's" not in file.parts}
            injected = bool(list(base.rglob("INJECTED")))
            return proc, native_args, runtime_cpu, injected, files

    @staticmethod
    def native_fixture(native_status, output_problem):
        # The fixture reads staged paths rather than treating argv as metadata.
        # S-mode emits all three adjusted files to exercise output collection.
        return "#!/usr/bin/env python3\n" + "\n".join([
            "import csv, json, os, pathlib, sys",
            "args = dict(zip(sys.argv[1::2], sys.argv[2::2]))",
            "root = pathlib.Path(__file__).resolve().parent.parent",
            "(root / 'native_args.json').write_text(json.dumps({'args': args, 'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS')}))",
            "print('>[Options] username: ' + args['--username'], flush=True)",
            "print('>[Options] token: ' + args['--token'], flush=True)",
            "print('Native fraction calculation started', flush=True)",
            "data = pathlib.Path(__file__).resolve().parent / 'data'",
            "for flag in ('--mixture', '--sigmatrix', '--refsample', '--sourceGEPs'):",
            "    if flag in args:",
            "        assert (data / args[flag]).is_file(), flag + ' is not staged'",
            "        assert (data / args[flag]).read_text(), flag + ' is empty'",
            "sys.exit(%d)" % native_status if native_status else "",
            "out = pathlib.Path(args['--outdir']).resolve()",
            "out.mkdir(exist_ok=True)",
            "samples = next(csv.reader((data / args['--mixture']).open(), delimiter='\\t'))[1:]",
            "rows = [['Mixture', 'B', 'CD4_T', 'P-value', 'Correlation', 'RMSE']] + [[sample, '0.25', '0.75', '0.01', '0.9', '0.1'] for sample in samples]",
            "problem = " + repr(output_problem),
            "if problem == 'empty': rows = []",
            "if problem == 'sample_order': rows = rows[:1] + rows[:0:-1]",
            "if problem == 'wrong_sample': rows[-1][0] = 'WRONG'",
            "if problem == 'wrong_cell': rows[0][1] = 'WRONG'",
            "if problem == 'nonnumeric': rows[1][1] = 'BAD'",
            "if problem == 'nonfinite': rows[1][1] = 'NaN'",
            "if problem == 'negative': rows[1][1:3] = ['-0.1', '1.1']",
            "if problem == 'row_sum': rows[1][1:3] = ['0.1', '0.1']",
            "name = 'CIBERSORTx_Adjusted.txt' if args.get('--rmbatchSmode') == 'TRUE' else 'CIBERSORTx_Results.txt'",
            "with (out / name).open('w', newline='') as stream:",
            "    csv.writer(stream, delimiter='\\t').writerows(rows)",
            "if args.get('--rmbatchSmode') == 'TRUE':",
            "    adjusted_mixture = (data / args['--mixture']).read_text()",
            "    if problem == 'adjusted_sample': adjusted_mixture = adjusted_mixture.replace('s3', 'WRONG', 1)",
            "    (out / 'CIBERSORTx_Mixtures_Adjusted.txt').write_text(adjusted_mixture)",
            "    (out / 'CIBERSORTx_sigmatrix_Adjusted.txt').write_text((data / args['--sigmatrix']).read_text() + 'not_in_mixture\\t3\\t4\\n')",
            "print('Native outputs written', flush=True)",
            "",
        ])

    def test_localized_inputs_and_literal_flags_reach_native_process(self):
        proc, native, cpu, injected, files = self.run_task()
        self.assertEqual(proc.returncode, 0, proc.stdout)
        args = native["args"]
        self.assertEqual(args["--username"], "user'$(touch INJECTED)@example.org")
        self.assertEqual(args["--token"], "token'`touch INJECTED`")
        self.assertEqual((cpu, native["OMP_NUM_THREADS"]), (8, "8"))
        self.assertEqual((args["--mixture"], args["--sigmatrix"]), ("mixture.txt", "signature.txt"))
        self.assertEqual((args["--perm"], args["--QN"]), ("100", "FALSE"))
        self.assertNotIn("--threads", args)  # This native Fractions build has no threads flag.
        self.assertNotIn("--refsample", args)
        self.assertNotIn("--sourceGEPs", args)
        self.assertFalse(injected)
        self.assertNotIn(args["--username"], proc.stdout)
        self.assertNotIn(args["--token"], proc.stdout)
        for field in ("stage=CIBERSORTxFractions", "start_time=", "completion_time=", "dimensions=", "outputs="):
            self.assertIn(field, proc.stdout)
        fractions = list(csv.reader(io.StringIO(files["results/fractions_only.txt"]), delimiter="\t"))
        self.assertEqual(fractions[0][1:], CELLS)
        self.assertEqual([row[0] for row in fractions[1:]], SAMPLES)
        self.assertTrue(all(len(row) == 3 for row in fractions))
        self.assertIn("results/signature_shared_genes.txt", files)

    def test_smode_stages_optional_files_and_returns_shared_gene_signature(self):
        proc, native, _, _, files = self.run_task(smode=True)
        self.assertEqual(proc.returncode, 0, proc.stdout)
        args = native["args"]
        self.assertEqual(args["--rmbatchSmode"], "TRUE")
        self.assertEqual((args["--refsample"], args["--sourceGEPs"]), ("refsample.txt", "source_geps.txt"))
        for name in ("CIBERSORTx_Adjusted.txt", "CIBERSORTx_Mixtures_Adjusted.txt", "CIBERSORTx_sigmatrix_Adjusted.txt"):
            self.assertIn("results/" + name, files)
        shared = files["results/signature_shared_genes.txt"]
        self.assertIn("sig\t1\t2", shared)
        self.assertNotIn("not_in_mixture", shared)

    def test_smode_requires_both_reference_files_before_native_process(self):
        for absent in ("refsample", "source_geps"):
            with self.subTest(absent=absent):
                proc, native, _, _, _ = self.run_task(lambda inputs, _: inputs.pop(absent), smode=True)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("Input error", proc.stdout)
                self.assertIsNone(native)

    def test_unresolved_cloud_files_fail_before_native_process(self):
        for name in ("mixture", "signature", "refsample", "source_geps"):
            with self.subTest(name=name):
                proc, native, _, _, _ = self.run_task(
                    lambda inputs, _, name=name: inputs.update({name: "gs://bucket/" + name + ".txt"}), smode=True)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("Localization error", proc.stdout)
                self.assertIsNone(native)

    def test_multiline_secret_is_literal_and_not_logged(self):
        token = "SECRET_FIRST'line\nFRACTIONS_TASK\nSECRET_LAST`touch INJECTED`"
        proc, native, _, injected, files = self.run_task(lambda inputs, _: inputs.update(token=token))
        self.assertEqual(proc.returncode, 0, proc.stdout)
        self.assertEqual(native["args"]["--token"], token)
        self.assertFalse(injected)
        for value in ("SECRET_FIRST", "FRACTIONS_TASK", "SECRET_LAST"):
            self.assertNotIn(value, proc.stdout)
            for name, content in files.items():
                if name.endswith(".log"):
                    self.assertNotIn(value, content)

    def test_native_nonzero_exit_is_not_hidden_by_logging(self):
        proc, _, _, _, _ = self.run_task(native_status=17)
        self.assertEqual(proc.returncode, 17, proc.stdout)
        self.assertIn("Native fraction calculation started", proc.stdout)

    def test_bad_fraction_outputs_fail_task(self):
        for problem in ("empty", "wrong_sample", "wrong_cell", "nonnumeric", "nonfinite", "negative", "row_sum"):
            with self.subTest(problem=problem):
                proc, _, _, _, _ = self.run_task(output_problem=problem)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("Output error", proc.stdout)

    def test_native_fraction_rows_are_reordered_to_input_samples(self):
        proc, _, _, _, files = self.run_task(output_problem="sample_order")
        self.assertEqual(proc.returncode, 0, proc.stdout)
        rows = list(csv.reader(io.StringIO(files["results/fractions_only.txt"]), delimiter="\t"))
        self.assertEqual([row[0] for row in rows[1:]], SAMPLES)

    def test_bad_input_numeric_value_fails_before_native_process(self):
        def invalid(inputs, _):
            Path(inputs["mixture"]).write_text("GeneSymbol\ts1\ts2\nsig\tNaN\t2\n")
        proc, native, _, _, _ = self.run_task(invalid)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("Input error", proc.stdout)
        self.assertIsNone(native)

    def test_no_gene_overlap_fails_before_native_process(self):
        # Dimensions alone do not establish overlap. Every reference that
        # participates in an S-mode fit must share genes with the mixture.
        for name in ("signature", "refsample", "source_geps"):
            with self.subTest(name=name):
                def no_overlap(inputs, _, name=name):
                    path = Path(inputs[name])
                    path.write_text(path.read_text().replace("sig\t", "absent_one\t").replace("ZNF804A\t", "absent_two\t"))
                proc, native, _, _, _ = self.run_task(no_overlap, smode=True)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("share no genes", proc.stdout)
                self.assertIsNone(native)

    def test_plain_shared_signature_excludes_genes_absent_from_mixture(self):
        def extra_signature_gene(inputs, _):
            path = Path(inputs["signature"])
            path.write_text(path.read_text() + "not_in_mixture\t3\t4\n")
        proc, _, _, _, files = self.run_task(extra_signature_gene)
        self.assertEqual(proc.returncode, 0, proc.stdout)
        rows = list(csv.reader(io.StringIO(files["results/signature_shared_genes.txt"]), delimiter="\t"))
        self.assertEqual(rows[0][1:], CELLS)
        self.assertEqual([row[0] for row in rows[1:]], ["sig"])

    def test_changed_adjusted_mixture_samples_fail_smode_task(self):
        proc, native, _, _, _ = self.run_task(output_problem="adjusted_sample", smode=True)
        self.assertIsNotNone(native)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("adjusted sample or cell labels changed", proc.stdout)

    def test_negative_permutation_count_fails_before_native_process(self):
        proc, native, _, _, _ = self.run_task(lambda inputs, _: inputs.update(permutations=-1))
        self.assertNotEqual(proc.returncode, 0)
        self.assertIsNone(native)

    def test_workflow_file_writer_check(self):
        path = ROOT / "scripts/check_wdl_file_scope.py"
        self.assertTrue(path.is_file(), "The workflow-scope regression check has not been created")
        spec = importlib.util.spec_from_file_location("check_fractions_workflow", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        forbidden = WDL.parse_document('version 1.0\nworkflow bad { Array[String] x = ["x"] File f = write_lines(x) }')
        self.assertTrue(module.workflow_file_writes(forbidden.workflow))
        nested = WDL.parse_document('version 1.0\ntask consumer { input { String text } command { echo ~{text} } }\nworkflow bad { call consumer { input: text=read_string(write_lines(["x"])) } }')
        self.assertTrue(module.workflow_file_writes(nested.workflow))
        workflow_output = WDL.parse_document('version 1.0\nworkflow bad { output { File text = write_lines(["x"]) } }')
        self.assertTrue(module.workflow_file_writes(workflow_output.workflow))
        scatter_decl = WDL.parse_document('version 1.0\nworkflow bad { scatter (x in ["x"]) { File text = write_lines([x]) } }')
        self.assertTrue(module.workflow_file_writes(scatter_decl.workflow))
        scatter_call = WDL.parse_document('version 1.0\ntask consumer { input { String text } command { echo ~{text} } }\nworkflow bad { scatter (x in ["x"]) { call consumer { input: text=read_string(write_lines([x])) } } }')
        self.assertTrue(module.workflow_file_writes(scatter_call.workflow))
        conditional_decl = WDL.parse_document('version 1.0\nworkflow bad { if (true) { File text = write_lines(["x"]) } }')
        self.assertTrue(module.workflow_file_writes(conditional_decl.workflow))
        task_only = WDL.parse_document('version 1.0\ntask ok { Array[String] x = ["x"] File f = write_lines(x) command { cat ~{f} } }\nworkflow good { call ok }')
        self.assertEqual(module.workflow_file_writes(task_only.workflow), [])
        self.assertEqual(module.workflow_file_writes(self.document().workflow), [])


if __name__ == "__main__":
    unittest.main()
