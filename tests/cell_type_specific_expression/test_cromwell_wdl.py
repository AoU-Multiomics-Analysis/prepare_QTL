"""Parse both public cell-type workflows with Cromwell's WDL validator.

MiniWDL accepts some syntax that Cromwell rejects. CI supplies a pinned
WOMTOOL_JAR; local runs may use the same jar and an optional WOMTOOL_JAVA path.
This check validates imports and types. It does not execute a Terra workflow.
"""
import os
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(os.environ.get("WOMTOOL_JAR"), "Set WOMTOOL_JAR to run Cromwell validation")
class CromwellWdlTest(unittest.TestCase):
    def test_public_workflows_and_imports_pass_cromwell_validation(self):
        jar = Path(os.environ["WOMTOOL_JAR"]).resolve()
        self.assertTrue(jar.is_file(), f"Womtool jar does not exist: {jar}")
        java = os.environ.get("WOMTOOL_JAVA", "java")
        for filename in ("deconvolution.wdl", "prepare_cell_type_eQTL.wdl"):
            with self.subTest(workflow=filename):
                result = subprocess.run(
                    [java, "-jar", str(jar), "validate",
                     str(ROOT / "workflows/cell_type_specific_expression" / filename)],
                    cwd=ROOT, capture_output=True, text=True, timeout=120,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
