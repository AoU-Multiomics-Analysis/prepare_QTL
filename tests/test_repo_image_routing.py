"""Validate declared, immutable routing for repository-owned WDL images."""
from pathlib import Path
import os
import re
import sys
import unittest

import WDL
import yaml


POLICY_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('RELEASE_SOURCE_ROOT', POLICY_ROOT)).resolve()
sys.path.insert(0, str(POLICY_ROOT / 'ci'))
from wdl_stage_routing import task_stages, validate_routing

WORKFLOWS = ROOT / "workflows"
DIGEST_REFERENCE = re.compile(r"^ghcr\.io/[^\s@]+@sha256:[0-9a-f]{64}$")

def literal_value(declaration):
    if declaration.expr is None or not isinstance(declaration.expr.literal, WDL.Value.String):
        return None
    return declaration.expr.literal.value


class RepositoryImageRoutingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load((POLICY_ROOT / "ci/image-stages.yml").read_text())
        cls.pins = yaml.safe_load((POLICY_ROOT / "ci/release-pins.yml").read_text())
        cls.repositories = {
            image: properties["repository"] for image, properties in cls.config["images"].items()
        }
        cls.documents = {}
        for path in sorted(WORKFLOWS.rglob("*.wdl")):
            relative = path.relative_to(ROOT).as_posix()
            if relative.startswith("workflows/genotype/"):
                continue
            cls.documents[relative] = WDL.load(str(path))

    def image_for_reference(self, value):
        matches = [image for image, repository in self.repositories.items()
                   if value.startswith(repository + "@")]
        self.assertEqual(len(matches), 1, value)
        return matches[0]

    def test_every_repository_runtime_uses_a_declared_string_input(self):
        for path, document in self.documents.items():
            for task in document.tasks:
                key = (path, task.name)
                stages = task_stages(task, self.config, ROOT)
                self.assertTrue(stages, key)
                self.assertIn("docker", task.runtime, key)
                image_input = str(task.runtime["docker"])
                declarations = {declaration.name: declaration for declaration in task.inputs}
                self.assertIn(image_input, declarations, key)
                self.assertEqual(str(declarations[image_input].type), "String", key)
                default = literal_value(declarations[image_input])
                if default is not None:
                    self.assertRegex(default, DIGEST_REFERENCE, key)
                    allowed_images = {self.config["stages"][stage]["image"] for stage in stages}
                    self.assertIn(self.image_for_reference(default), allowed_images, key)

    def test_calls_forward_the_matching_workflow_image(self):
        self.assertEqual(validate_routing(ROOT, POLICY_ROOT), [])

    def test_every_image_default_is_registered_and_immutable(self):
        configured = set()
        for stage, targets in self.pins["stages"].items():
            repository = self.repositories[self.config["stages"][stage]["image"]]
            for target in targets:
                scope = target.get("scope", "workflow")
                key = (target["path"], scope, target.get("task"), target["input"])
                self.assertNotIn(key, configured)
                configured.add(key)
                document = WDL.load(str(ROOT / target["path"]))
                if scope == "workflow":
                    self.assertIsNotNone(document.workflow, key)
                    declarations = document.workflow.inputs
                else:
                    self.assertEqual(scope, "task", key)
                    matching_tasks = [task for task in document.tasks if task.name == target.get("task")]
                    self.assertEqual(len(matching_tasks), 1, key)
                    declarations = matching_tasks[0].inputs
                matching = [declaration for declaration in declarations
                            if declaration.name == target["input"]]
                self.assertEqual(len(matching), 1, key)
                value = literal_value(matching[0])
                self.assertIsNotNone(value, key)
                self.assertRegex(value, DIGEST_REFERENCE, key)
                self.assertTrue(value.startswith(repository + "@"), key)

        actual_defaults = set()
        for path, document in self.documents.items():
            if document.workflow is not None:
                for declaration in document.workflow.inputs:
                    value = literal_value(declaration)
                    if value and value.startswith("ghcr.io/aou-multiomics-analysis/"):
                        actual_defaults.add((path, "workflow", None, declaration.name))
            for task in document.tasks:
                for declaration in task.inputs:
                    value = literal_value(declaration)
                    if value and value.startswith("ghcr.io/aou-multiomics-analysis/"):
                        actual_defaults.add((path, "task", task.name, declaration.name))
        self.assertTrue(actual_defaults.issubset(configured), sorted(actual_defaults - configured))


if __name__ == "__main__":
    unittest.main()
