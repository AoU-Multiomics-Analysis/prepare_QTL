"""Verify the canonical model dependency at each public workflow boundary."""
import unittest
from pathlib import Path
import WDL

ROOT = Path(__file__).resolve().parents[2]


class TcaCleanupTest(unittest.TestCase):
    def test_consumers_and_public_output_use_cleaned_model(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/deconvolution.wdl"))
        calls = {n.name: n for n in doc.workflow.body if isinstance(n, WDL.Tree.Call)}
        self.assertIn("CleanTcaModel", calls)
        self.assertEqual(str(calls["CleanTcaModel"].inputs["unfiltered_model"]),
                         "select_first([precomputed_tca_model, FitTca.model])")
        for consumer in ("ExportTcaBeds", "BuildManifest"):
            self.assertEqual(str(calls[consumer].inputs["model"]), "CleanTcaModel.model")
        outputs = {n.name: n for n in doc.workflow.outputs}
        for name, expression in {
            "tca_model": "CleanTcaModel.model",
            "tca_model_unfiltered": "FitTca.model",
            "tca_numerical_excluded_genes": "CleanTcaModel.excluded_genes",
            "tca_cleanup_log": "CleanTcaModel.log",
        }.items():
            self.assertEqual(str(outputs[name].expr), expression)

    def test_eqtl_entry_point_exposes_final_and_audit_outputs(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"))
        outputs = {n.name: n for n in doc.workflow.outputs}
        for name in ("tca_model", "tca_model_unfiltered", "tca_numerical_excluded_genes", "tca_cleanup_log"):
            self.assertIn(name, outputs)
            self.assertEqual(str(outputs[name].expr), "CellTypeDeconvolution." + name)


if __name__ == "__main__":
    unittest.main()
