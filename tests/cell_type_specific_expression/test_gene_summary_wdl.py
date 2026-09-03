"""Verify that summaries consume exported CPM, not QTL-transformed phenotypes."""
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]


class GeneSummaryTest(unittest.TestCase):
    def test_summary_consumes_exported_files_and_is_a_public_output(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/deconvolution.wdl"))
        calls = {n.name: n for n in doc.workflow.body if isinstance(n, WDL.Tree.Call)}
        self.assertIn("SummarizeCellTypeBeds", calls)
        summary = calls["SummarizeCellTypeBeds"]
        self.assertEqual(str(summary.inputs["cell_type_beds"]), "ExportTcaBeds.cell_type_beds")
        self.assertEqual(str(summary.inputs["cell_type_bed_inventory"]),
                         "ExportTcaBeds.cell_type_bed_inventory")
        outputs = {n.name: n for n in doc.workflow.outputs}
        self.assertEqual(str(outputs["cell_type_gene_summary"].expr),
                         "SummarizeCellTypeBeds.summary")
        self.assertEqual(str(outputs["gene_summary_log"].type), "File")

    def test_eqtl_entry_point_exposes_same_summary(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"))
        outputs = {n.name: n for n in doc.workflow.outputs}
        self.assertIn("cell_type_gene_summary", outputs)
        self.assertEqual(str(outputs["cell_type_gene_summary"].expr),
                         "CellTypeDeconvolution.cell_type_gene_summary")


if __name__ == "__main__":
    unittest.main()
