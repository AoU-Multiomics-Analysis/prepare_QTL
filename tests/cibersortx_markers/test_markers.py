"""Test the real WDL command with a stand-in for the licensed calculation.

File localization and output upload are simulated. These tests do not run Terra
or establish the accuracy of the CIBERSORTx marker selection algorithm. Set
CIBERSORTX_TEST_IMAGE in CI to execute the same command in the existing image.
"""
import glob
import hashlib
import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]
WDL_PATH = ROOT / "workflows/cell_type_specific_expression/cibersortx_markers.wdl"
REFERENCE = (
    "GeneSymbol\tB\tB\tB\tB\tB\tCD4_T\tCD4_T\tCD4_T\tCD4_T\tCD4_T\n"
    "geneA\t10\t11\t12\t13\t14\t1\t1\t1\t1\t1\n"
    "geneB\t1\t1\t1\t1\t1\t10\t11\t12\t13\t14\n"
)
USERNAME = "user'$(touch INJECTED)@example.org"
TOKEN = "token'`touch INJECTED`"


class OutputFunctions(WDL.StdLib.TaskOutputs):
    """Resolve task output files from the temporary execution directory."""

    def __init__(self, base):
        super().__init__("1.0")
        self.stdout = WDL.StdLib.StaticFunction(
            "stdout", [], WDL.Type.File(), lambda: WDL.Value.File(str(base / "stdout.txt")))
        self.stderr = WDL.StdLib.StaticFunction(
            "stderr", [], WDL.Type.File(), lambda: WDL.Value.File(str(base / "stderr.txt")))
        self.glob = WDL.StdLib.StaticFunction(
            "glob", [WDL.Type.String()], WDL.Type.Array(WDL.Type.File()),
            lambda pattern: WDL.Value.Array(WDL.Type.File(), [
                WDL.Value.File(path) for path in sorted(glob.glob(str(base / pattern.value)))
            ]))


class MarkerTest(unittest.TestCase):
    def document(self):
        self.assertTrue(WDL_PATH.is_file(), "The marker derivation WDL has not been created")
        return WDL.load(str(WDL_PATH))

    def run_task(self, change=None, native_status=0, output_mode="valid"):
        doc = self.document()
        task = next(task for task in doc.tasks if task.name == "DeriveMarkers")
        with tempfile.TemporaryDirectory(prefix="markers localized ") as tmp:
            base = Path(tmp)
            localized = base / "localized input's"
            localized.mkdir()
            reference = localized / "reference ' $(touch INJECTED).tsv"
            reference.write_text(REFERENCE)
            inputs = {"reference": str(reference), "username": USERNAME, "token": TOKEN}
            if change:
                change(inputs, localized)

            # The engine sees cloud Files before localization. String inputs
            # must stay untouched when these File values become local paths.
            cloud_inputs = dict(inputs)
            downloads = {}
            if "://" not in inputs["reference"]:
                uri = "gs://test-input-bucket/reference.tsv"
                downloads[uri] = inputs["reference"]
                cloud_inputs["reference"] = uri
            self.assertIsInstance(doc.workflow.available_inputs.resolve("reference").type, WDL.Type.File)
            self.assertIsInstance(task.available_inputs.resolve("reference").type, WDL.Type.File)
            env = WDL.values_from_json(cloud_inputs, task.available_inputs)
            stdlib = WDL.StdLib.Base("1.0")
            for decl in task.inputs:
                if not env.has_binding(decl.name) and decl.expr is not None:
                    env = env.bind(decl.name, decl.expr.eval(env, stdlib))
            # Evaluate task declarations before input localization, as Terra
            # can do. A declaration that serializes a cloud path is unsafe.
            for decl in task.postinputs:
                env = env.bind(decl.name, decl.expr.eval(env, stdlib))
            env = WDL.Value.rewrite_env_paths(env, lambda file: downloads.get(file.value, file.value))
            runtime_cpu = task.runtime["cpu"].eval(env, stdlib).value
            command = task.command.eval(env, stdlib).value

            src = base / "src"
            src.mkdir()
            native = src / "CIBERSORTxFractions"
            image = os.environ.get("CIBERSORTX_TEST_IMAGE")
            data_path = "/src/data" if image else str(src / "data")
            native.write_text(self.native_fixture(base, data_path, native_status, output_mode))
            native.chmod(0o755)
            if not image:
                command = command.replace("/src", str(src))
            script = base / "command.sh"
            script.write_text(command)
            invocation = ["bash", str(script)]
            if image:
                invocation = [
                    "docker", "run", "--rm", "--platform", "linux/amd64",
                    "--entrypoint", "/bin/bash", "--workdir", str(base),
                    "-v", str(base) + ":" + str(base),
                    "-v", str(native) + ":/src/CIBERSORTxFractions:ro", image,
                    "-c", 'set +e; bash "$1"; result=$?; chown -R "$2:$3" "$4"; exit "$result"',
                    "markers-test", str(script), str(os.getuid()), str(os.getgid()), str(base),
                ]
            proc = subprocess.run(invocation, cwd=base, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=45)
            (base / "stdout.txt").write_text(proc.stdout)
            (base / "stderr.txt").write_text(proc.stderr)
            args_path = base / "native_args.json"
            report_paths = [base / name for name in ("input_validation.json", "output_validation.json")]
            reports = [json.loads(path.read_text()) if path.exists() else None for path in report_paths]
            outputs = {}
            workflow_outputs = {}
            if proc.returncode == 0:
                output_env = env
                output_stdlib = OutputFunctions(base)
                for decl in task.outputs:
                    value = decl.expr.eval(output_env, output_stdlib).coerce(decl.type)
                    value = WDL.Value.rewrite_paths(value, lambda file: str(base / file.value))
                    output_env = output_env.bind(decl.name, value)
                    outputs[decl.name] = value
                # Upload completed task Files, then pass their cloud File
                # values through the actual workflow output expressions.
                flow_env = WDL.Env.Bindings()
                for name, value in outputs.items():
                    uploaded = WDL.Value.rewrite_paths(
                        value, lambda file: "gs://test-output-bucket/" + Path(file.value).name)
                    flow_env = flow_env.bind("DeriveMarkers." + name, uploaded)
                for decl in doc.workflow.outputs:
                    workflow_outputs[decl.name] = decl.expr.eval(flow_env, stdlib).coerce(decl.type)
            result = {
                "process": proc,
                "args": json.loads(args_path.read_text()) if args_path.exists() else None,
                "cpu": runtime_cpu,
                "injected": bool(list(base.rglob("INJECTED"))),
                "input_report": reports[0], "output_report": reports[1],
                "workflow_outputs": workflow_outputs,
                "output_files": {name: Path(value.value).read_text() for name, value in outputs.items()
                                 if isinstance(value, WDL.Value.File) and Path(value.value).is_file()},
                "run_log": (base / "cibersortx_markers.log").read_text()
                           if (base / "cibersortx_markers.log").exists() else "",
            }
            return result

    @staticmethod
    def native_fixture(base, data_path, status, mode):
        # Only the authenticated native calculation is replaced. The WDL
        # still validates inputs, stages files, logs, and validates outputs.
        return "#!/usr/bin/env python3\n" + f"""
import hashlib, json, pathlib, sys
args = dict(zip(sys.argv[1::2], sys.argv[2::2]))
out = pathlib.Path(args['--outdir']).resolve()
reference = pathlib.Path({data_path!r}) / args['--refsample']
args['_reference_sha256'] = hashlib.sha256(reference.read_bytes()).hexdigest()
pathlib.Path({str(base / 'native_args.json')!r}).write_text(json.dumps(args))
print('Native calculation started', flush=True)
print('>[Options] username: ' + args['--username'], flush=True)
print('>[Options] token: ' + args['--token'], flush=True)
print('Native credential check: ' + args['--token'], file=sys.stderr, flush=True)
if {status!r}:
    sys.exit({status!r})
mode = {mode!r}
if mode == 'no_output':
    sys.exit(0)
replicates = int(args['--replicates'])
types = ['B', 'CD4_T']
profile_headers = [cell + ('.' + str(i) if i else '') for cell in types for i in range(replicates)]
signature = 'NAME\\tB\\tCD4_T\\ngeneA\\t10\\t1\\ngeneB\\t1\\t10\\n'
source = 'genesymbols\\tB\\tCD4_T\\ngeneA\\t10\\t1\\ngeneB\\t1\\t10\\n'
profiles = 'GeneSymbol\\t' + '\\t'.join(profile_headers) + '\\n'
profiles += 'geneA\\t' + '\\t'.join(['10'] * replicates + ['1'] * replicates) + '\\n'
profiles += 'geneB\\t' + '\\t'.join(['1'] * replicates + ['10'] * replicates) + '\\n'
classes = 'B\\t' + '\\t'.join(['1'] * replicates + ['2'] * replicates) + '\\n'
classes += 'CD4_T\\t' + '\\t'.join(['2'] * replicates + ['1'] * replicates) + '\\n'
if mode == 'negative_signature':
    signature = signature.replace('geneA\\t10', 'geneA\\t-1')
elif mode == 'nonfinite_source':
    source = source.replace('geneA\\t10', 'geneA\\tnan')
elif mode == 'signature_cell_order':
    signature = signature.replace('NAME\\tB\\tCD4_T', 'NAME\\tCD4_T\\tB')
elif mode == 'signature_gene_absent':
    signature = signature.replace('geneA', 'absent_gene')
elif mode == 'truncated_reference':
    profiles = '\\n'.join(line.rsplit('\\t', 1)[0] for line in profiles.splitlines()) + '\\n'
elif mode == 'invalid_phenotype':
    classes = classes.replace('B\\t1', 'B\\t0', 1)
prefix = 'CIBERSORTx_reference_inferred_phenoclasses.CIBERSORTx_reference_inferred_refsample.bm.K'
(out / (prefix + args['--k.max'] + '.txt')).write_text(signature)
if mode != 'missing_source':
    (out / 'CIBERSORTx_cell_type_sourceGEP.txt').write_text(source)
(out / 'CIBERSORTx_reference_inferred_refsample.txt').write_text(profiles)
(out / 'CIBERSORTx_reference_inferred_phenoclasses.txt').write_text(classes)
if mode == 'pdf':
    (out / (prefix + args['--k.max'] + '.pdf')).write_bytes(b'%PDF-1.4\\nfixture\\n')
"""

    def assert_input_rejected(self, change):
        result = self.run_task(change)
        self.assertNotEqual(result["process"].returncode, 0)
        self.assertIsNone(result["args"], result["process"].stdout + result["process"].stderr)
        self.assertIn("error", (result["process"].stdout + result["process"].stderr).lower())

    def test_localized_reference_and_literal_credentials_reach_native(self):
        # Catch unsafe shell interpolation and failures to stage the File.
        result = self.run_task()
        proc, args = result["process"], result["args"]
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertEqual(args["--refsample"], "reference.tsv")
        self.assertEqual(args["_reference_sha256"], hashlib.sha256(REFERENCE.encode()).hexdigest())
        self.assertEqual(args["--username"], USERNAME)
        self.assertEqual(args["--token"], TOKEN)
        self.assertFalse(result["injected"])
        self.assertEqual(args["--single_cell"], "TRUE")
        self.assertEqual(args["--QN"], "FALSE")
        self.assertTrue({"--mixture", "--sigmatrix", "--rmbatchSmode", "--rmbatchBmode"}.isdisjoint(args))
        self.assertEqual(result["input_report"]["gene_count"], 2)
        self.assertEqual(result["input_report"]["cell_count"], 10)
        self.assertEqual(result["input_report"]["cell_types"], ["B", "CD4_T"])
        self.assertEqual(result["output_report"]["gene_count"], 2)
        self.assertEqual(result["output_report"]["reference_profile_count"], 10)

    def test_parameter_overrides_reach_native(self):
        # Catch ignored user settings or a CPU request that does not change.
        result = self.run_task(lambda inputs, _: inputs.update(
            replicates=3, sampling=0.5, fraction=0.2, min_genes=1,
            max_genes=2, q_value=0.02, max_condition_number=100, cpu=2))
        self.assertEqual(result["process"].returncode, 0,
                         result["process"].stdout + result["process"].stderr)
        args = result["args"]
        self.assertEqual(args["--replicates"], "3")
        self.assertEqual(float(args["--sampling"]), 0.5)
        self.assertEqual(float(args["--fraction"]), 0.2)
        self.assertEqual(args["--G.min"], "1")
        self.assertEqual(args["--G.max"], "2")
        self.assertEqual(float(args["--q.value"]), 0.02)
        self.assertEqual(args["--k.max"], "100")
        self.assertEqual(result["cpu"], 2)
        self.assertEqual(result["output_report"]["reference_profile_count"], 6)

    def test_outputs_stay_files_after_cloud_upload(self):
        # Catch output paths becoming Strings before the workflow consumer.
        result = self.run_task(output_mode="pdf")
        self.assertEqual(result["process"].returncode, 0,
                         result["process"].stdout + result["process"].stderr)
        expected = {"signature_matrix": "signature_matrix.tsv", "source_geps": "source_geps.tsv",
                    "reference_sample": "reference_sample.tsv", "phenotype_classes": "phenotype_classes.tsv"}
        for name, filename in expected.items():
            with self.subTest(output=name):
                self.assertIn(name, result["output_files"])
                value = result["workflow_outputs"][name]
                self.assertIsInstance(value, WDL.Value.File)
                self.assertEqual(value.value, "gs://test-output-bucket/" + filename)

    def test_unresolved_cloud_file_fails_before_native(self):
        result = self.run_task(lambda inputs, _: inputs.update(reference="gs://bucket/reference.tsv"))
        self.assertNotEqual(result["process"].returncode, 0)
        self.assertIn("Localization error", result["process"].stdout + result["process"].stderr)
        self.assertIsNone(result["args"])

    def test_unreadable_file_fails_before_native(self):
        self.assert_input_rejected(lambda inputs, _: Path(inputs["reference"]).unlink())

    def test_malformed_references_fail_before_native(self):
        # Each case changes a property that the native parser relies on.
        invalid = {
            "negative": REFERENCE.replace("geneA\t10", "geneA\t-1"),
            "nonfinite": REFERENCE.replace("geneA\t10", "geneA\tNaN"),
            "nonnumeric": REFERENCE.replace("geneA\t10", "geneA\tbad"),
            "unequal_width": REFERENCE.replace("geneA\t10\t11", "geneA\t10"),
            "duplicate_gene": REFERENCE.replace("geneB", "geneA"),
            "one_cell_type": REFERENCE.replace("CD4_T", "B"),
            "too_few_cells": "GeneSymbol\tB\tB\tCD4_T\tCD4_T\ngeneA\t1\t2\t3\t4\n",
            "empty": "",
        }
        for name, contents in invalid.items():
            with self.subTest(case=name):
                self.assert_input_rejected(lambda inputs, _, text=contents: Path(inputs["reference"]).write_text(text))

    def test_invalid_parameters_fail_before_native(self):
        invalid = [
            {"replicates": 1}, {"sampling": 0.0}, {"sampling": 1.1},
            {"fraction": -0.1}, {"fraction": 1.1}, {"min_genes": 0},
            {"min_genes": 5, "max_genes": 4}, {"q_value": 0.0},
            {"q_value": 1.1}, {"max_condition_number": 0},
            {"cpu": 0}, {"memory_gb": 0}, {"disk_gb": 0},
        ]
        for changed in invalid:
            with self.subTest(parameters=changed):
                self.assert_input_rejected(lambda inputs, _, values=changed: inputs.update(values))

    def test_native_failure_survives_logging(self):
        result = self.run_task(native_status=17)
        self.assertNotEqual(result["process"].returncode, 0)
        self.assertIn("exit status 17", result["run_log"])
        self.assertIn("Native calculation started", result["run_log"])

    def test_missing_and_malformed_native_outputs_fail(self):
        # Native exit zero is insufficient when its files cannot be consumed.
        for mode in ("no_output", "missing_source", "negative_signature", "nonfinite_source",
                     "signature_cell_order", "signature_gene_absent", "truncated_reference", "invalid_phenotype"):
            with self.subTest(output=mode):
                result = self.run_task(output_mode=mode)
                self.assertIsNotNone(result["args"])
                self.assertNotEqual(result["process"].returncode, 0)
                self.assertIn("Output error", result["process"].stdout + result["process"].stderr)

    def test_credentials_are_removed_from_all_logs(self):
        result = self.run_task()
        self.assertEqual(result["process"].returncode, 0,
                         result["process"].stdout + result["process"].stderr)
        for log in (result["process"].stdout, result["process"].stderr, result["run_log"]):
            with self.subTest(log=log[:40]):
                self.assertNotIn(USERNAME, log)
                self.assertNotIn(TOKEN, log)
        self.assertIn("Native calculation started", result["run_log"])

    def test_multiline_credential_remains_literal_and_redacted(self):
        token = "first'line\nMARKERS_TASK\nlast`touch INJECTED`"
        result = self.run_task(lambda inputs, _: inputs.update(token=token))
        self.assertEqual(result["process"].returncode, 0,
                         result["process"].stdout + result["process"].stderr)
        self.assertEqual(result["args"]["--token"], token)
        self.assertFalse(result["injected"])
        for log in (result["process"].stdout, result["process"].stderr, result["run_log"]):
            for secret_line in token.splitlines():
                self.assertNotIn(secret_line, log)

    def test_workflow_has_no_engine_scope_file_writers(self):
        # Use the shared regression checker rather than a text search.
        spec = importlib.util.spec_from_file_location("check_scope", ROOT / "scripts/check_wdl_file_scope.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.workflow_file_writes(self.document().workflow), [])


if __name__ == "__main__":
    unittest.main()
