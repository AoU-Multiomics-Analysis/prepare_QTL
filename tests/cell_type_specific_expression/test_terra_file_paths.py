"""Regression checks for optional File localization and Terra workflow scope.

The command test models evaluation before localization, then rewrites only File
values. It is not an execution of Cromwell or a complete Terra workflow.
"""
import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]


def render_localized_command(task, optional_name, local_path, present):
    stdlib = WDL.StdLib.Base("1.0")
    bindings = WDL.Env.Bindings()
    for declaration in task.inputs:
        if declaration.name == optional_name:
            value = WDL.Value.File("gs://fixture/remote.json") if present else WDL.Value.Null()
        elif declaration.expr is not None:
            value = declaration.expr.eval(bindings, stdlib)
        elif isinstance(declaration.type, WDL.Type.File):
            value = WDL.Value.File("gs://fixture/input.tsv")
        elif isinstance(declaration.type, WDL.Type.Array):
            value = WDL.Value.Array(WDL.Type.String(), [WDL.Value.String("protein_coding")])
        elif isinstance(declaration.type, WDL.Type.Float):
            value = WDL.Value.Float(0.1)
        elif isinstance(declaration.type, WDL.Type.Int):
            value = WDL.Value.Int(1)
        elif isinstance(declaration.type, WDL.Type.Boolean):
            value = WDL.Value.Boolean(False)
        else:
            value = WDL.Value.String("fixture")
        bindings = bindings.bind(declaration.name, value.coerce(declaration.type))
    for declaration in task.postinputs:
        bindings = bindings.bind(declaration.name,
                                 declaration.expr.eval(bindings, stdlib).coerce(declaration.type))
    required_path = local_path.parent / "required.tsv"
    required_path.write_text("header\nrow\n")
    bindings = WDL.Value.rewrite_env_paths(bindings, lambda value: str(
        local_path if value.value == "gs://fixture/remote.json" else required_path))
    return task.command.eval(bindings, stdlib).value


class OptionalFileLocalizationTest(unittest.TestCase):
    def test_manifest_reports_unreadable_local_metadata_before_running_r(self):
        task = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/tasks/qc.wdl")).tasks[0]
        with tempfile.TemporaryDirectory() as directory:
            command = render_localized_command(task, "hspe_metadata",
                                                Path(directory) / "missing.json", True)
            result = subprocess.run(["bash"], input=command, text=True,
                                    cwd=directory, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("hspe_metadata_not_readable", result.stderr)

    def test_commands_pass_localized_optional_files_or_omit_the_argument(self):
        cases = [("qc.wdl", "BuildManifest", "hspe_metadata", "--hspe-metadata"),
                 ("tca.wdl", "FitTca", "covariates", "--covariates"),
                 ("tca.wdl", "ExportTcaBeds", "covariates", "--covariates")]
        for source, task_name, name, flag in cases:
            task = next(t for t in WDL.load(str(ROOT / "workflows/cell_type_specific_expression/tasks" / source)).tasks
                        if t.name == task_name)
            for present in (True, False):
                with self.subTest(task=task_name, present=present), tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    local_path = root / "donor's localized metadata $literal `echo unsafe` $(echo unsafe).json"
                    local_path.write_text('{"source":"localized"}\n')
                    command = render_localized_command(task, name, local_path, present)
                    # Replace only the expensive R program; inspect the actual shell argv.
                    # Do not replace any command interpolation or optional-argument logic.
                    capture = '''
Rscript() { python3 -c 'import json,sys; print(json.dumps(sys.argv[1:]))' "$@" > argv.json; }
export -f Rscript
mkdir -p outputs
printf 'header\\nrow\\n' > outputs/cell_type_bed_inventory.tsv
'''
                    result = subprocess.run(["bash"], input=capture + command, text=True,
                                            cwd=root, capture_output=True, env=os.environ.copy())
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    argv = json.loads((root / "argv.json").read_text())
                    if present:
                        self.assertIn(flag, argv)
                        self.assertEqual(argv[argv.index(flag) + 1], str(local_path))
                        self.assertEqual(json.loads(Path(argv[argv.index(flag) + 1]).read_text()),
                                         {"source": "localized"})
                    else:
                        self.assertNotIn(flag, argv)
                    self.assertNotIn("gs://", command)


class WorkflowScopeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / "scripts/check_wdl_file_scope.py"
        if not path.exists():
            raise AssertionError("The Terra workflow file-scope checker is missing")
        spec = importlib.util.spec_from_file_location("file_scope", path)
        cls.checker = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.checker)

    def test_checks_inputs_calls_outputs_and_nested_workflow_expressions(self):
        cases = [
            'input { File f = write_lines(["sample"]) }',
            'File f = write_lines(["sample"])',
            'call sink { input: f = write_tsv([["sample"]]) }',
            'output { File f = write_map({"sample": "path"}) }',
            'scatter (s in ["sample"]) { if (true) { call sink { input: f = write_json([s]) } } }',
            'scatter (s in read_lines(write_lines(["sample"]))) { String x = s }',
            'if (size(write_lines(["sample"])) > 0) { String x = "sample" }',
        ]
        for body in cases:
            with self.subTest(body=body):
                document = WDL.parse_document('version 1.0\n'
                    'task sink { input { File f } command { echo "~{f}" } }\n'
                    f'workflow example {{ {body} }}')
                self.assertEqual(len(self.checker.workflow_file_writes(document.workflow)), 1)

    def test_task_file_writes_are_allowed_and_comments_are_not_expressions(self):
        document = WDL.parse_document('version 1.0\n'
            'task sink { input { String s } File f = write_lines([s]) '
            'command { cat "~{f}" "~{write_json([s])}" } }\n'
            'workflow example { # write_lines(["not executed"])\n'
            'call sink { input: s = "sample" } }')
        self.assertEqual(self.checker.workflow_file_writes(document.workflow), [])

    def test_repo_workflows_and_checker_exit_status(self):
        result = subprocess.run([os.sys.executable, str(ROOT / "scripts/check_wdl_file_scope.py"),
                                 str(ROOT / "workflows")], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with tempfile.TemporaryDirectory() as directory:
            invalid = Path(directory) / "invalid.wdl"
            invalid.write_text('version 1.0\nworkflow w { File f = write_json(["x"]) }')
            result = subprocess.run([os.sys.executable, str(ROOT / "scripts/check_wdl_file_scope.py"),
                                     str(invalid)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("write_json", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
