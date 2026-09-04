"""Check expression scale routing and optional-file command arguments."""
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[1]
STDLIB = WDL.StdLib.Base("1.0")


def calls_in(nodes):
    for node in nodes:
        if isinstance(node, WDL.Tree.Call):
            yield node
        elif isinstance(node, (WDL.Tree.Scatter, WDL.Tree.Conditional)):
            yield from calls_in(node.body)


class ExpressionInputTest(unittest.TestCase):
    def test_tca_beds_reach_the_linear_cpm_input(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"))
        call = next(c for c in calls_in(doc.workflow.body) if c.name == "PrepareCellTypeEqtl")
        self.assertIn("CpmBed", call.inputs)
        self.assertNotIn("Log2CpmBed", call.inputs)
        bindings = WDL.Env.Bindings().bind("index", WDL.Value.Int(1)).bind(
            "PrepareScatterInputs.expression_beds",
            WDL.Value.Array(WDL.Type.File(), [WDL.Value.File("b.bed.gz"), WDL.Value.File("cd4.bed.gz")]))
        self.assertEqual(call.inputs["CpmBed"].eval(bindings, STDLIB).value, "cd4.bed.gz")

    def test_all_expression_modes_reach_the_r_command_without_changing_paths(self):
        doc = WDL.load(str(ROOT / "workflows/expression/prepare_eQTL.wdl"))
        call = next(c for c in calls_in(doc.workflow.body) if c.name == "eqtl_prepare_expression")
        task = call.callee
        modes = ("CountGCT", "CpmBed", "Log2CpmBed")
        # This test receives already-localized File values, as command rendering does.
        for mode in modes:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                path = str(Path(directory) / ("donor's $literal `echo unsafe`.bed.gz" if mode == "CpmBed" else "input bed.gz"))
                bindings = WDL.Env.Bindings()
                for declaration in task.inputs:
                    name = declaration.name
                    if name in modes:
                        value = WDL.Value.File(path) if name == mode else WDL.Value.Null()
                    elif name == "AnnotationGTF":
                        value = WDL.Value.File("genes.gtf") if mode == "CountGCT" else WDL.Value.Null()
                    elif declaration.expr is not None:
                        value = declaration.expr.eval(bindings, STDLIB)
                    elif isinstance(declaration.type, WDL.Type.File):
                        value = WDL.Value.File("samples.tsv")
                    elif isinstance(declaration.type, WDL.Type.Int):
                        value = WDL.Value.Int(1)
                    else:
                        value = WDL.Value.String("fixture")
                    bindings = bindings.bind(name, value.coerce(declaration.type))
                self.assertIn(mode, call.inputs)
                forwarded = call.inputs[mode].eval(bindings, STDLIB)
                self.assertEqual(forwarded.value, path)
                command = task.command.eval(bindings, STDLIB).value
                capture = "Rscript() { python3 -c 'import json,sys; print(json.dumps(sys.argv[1:]))' \"$@\" > argv.json; }\n"
                result = subprocess.run(["bash"], input=capture + command, text=True,
                                        cwd=directory, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                argv = json.loads((Path(directory) / "argv.json").read_text())
                self.assertEqual(argv[argv.index("--" + mode) + 1], path)
                for other in set(modes) - {mode}:
                    self.assertNotIn("--" + other, argv)


if __name__ == "__main__":
    unittest.main()
