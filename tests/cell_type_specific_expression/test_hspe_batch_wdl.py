"""Check the typed WDL boundaries of the HSPE scatter."""
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]


class HspeScatterTest(unittest.TestCase):
    def test_workers_receive_small_files_and_merge_precedes_group_filtering(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/deconvolution.wdl"))
        condition = next(n for n in doc.workflow.body if isinstance(n, WDL.Tree.Conditional))
        calls = {n.name: n for n in condition.body if isinstance(n, WDL.Tree.Call)}
        self.assertIn("PrepareHspeBatches", calls)
        self.assertIn("MergeHspeBatches", calls)
        scatter = next(n for n in condition.body if isinstance(n, WDL.Tree.Scatter))
        worker = next(n for n in scatter.body if isinstance(n, WDL.Tree.Call))
        self.assertEqual(worker.name, "RunHspeBatch")
        worker_files = {n.name for n in worker.callee.inputs if isinstance(n.type, WDL.Type.File)}
        self.assertEqual(worker_files, {"prepared", "batch"})
        self.assertEqual(str(calls["MergeHspeBatches"].inputs["batch_results"].type), "Array[File]")
        selected = next(n for n in doc.workflow.body
                        if isinstance(n, WDL.Tree.Decl) and n.name == "proportions_for_processing")
        self.assertIn("MergeHspeBatches.proportions", str(selected.expr))

    def test_batch_settings_are_exposed_and_forwarded(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"))
        call = next(n for n in doc.workflow.body
                    if isinstance(n, WDL.Tree.Call) and n.name == "CellTypeDeconvolution")
        for name in ("hspe_batch_size", "hspe_batch_memory", "hspe_batch_disk_gb"):
            self.assertIn(name, {n.name for n in doc.workflow.inputs})
            self.assertIn(name, call.inputs)
        outputs = {n.name: str(n.type) for n in doc.workflow.outputs}
        self.assertEqual(outputs["hspe_sample_diagnostics"], "File?")


if __name__ == "__main__":
    unittest.main()
