"""Validate declared, immutable routing for repository-owned WDL images."""
from pathlib import Path
import os
import re
import unittest

import WDL
import yaml


POLICY_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('RELEASE_SOURCE_ROOT', POLICY_ROOT)).resolve()
WORKFLOWS = ROOT / "workflows"
DIGEST_REFERENCE = re.compile(r"^ghcr\.io/[^\s@]+@sha256:[0-9a-f]{64}$")

# Each non-genotype runtime in these workflow families uses an image that this
# repository builds. This literal table makes a newly hard-coded or unclassified
# runtime fail until its release family is selected deliberately.
TASK_IMAGES = {
    ("workflows/cell_type_specific_expression/tasks/expression.wdl", "FilterExpressionGenes"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/gene_summary.wdl", "SummarizeCellTypeBeds"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/hspe.wdl", "PrepareHspeBatches"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/hspe.wdl", "RunHspeBatch"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/hspe.wdl", "MergeHspeBatches"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/integration.wdl", "PrepareScatterInputs"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/integration.wdl", "BuildQtlManifest"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/proportions.wdl", "ValidateProportionMode"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/proportions.wdl", "ProcessProportions"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/qc.wdl", "BuildManifest"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/reference_filter.wdl", "PrepareHaemopedia"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/reference_filter.wdl", "FilterCellTypeBeds"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/tca.wdl", "FitTca"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/tca.wdl", "CleanTcaModel"): "cell_type",
    ("workflows/cell_type_specific_expression/tasks/tca.wdl", "ExportTcaBeds"): "cell_type",
    ("workflows/common/MergeCovariates.wdl", "MergeCovariatesR"): "standard",
    ("workflows/common/ResidualizePhenotypes.wdl", "ResidualizePhenotypes"): "standard",
    ("workflows/common/calculate_phenotypePCs.wdl", "ComputePCs"): "standard",
    ("workflows/expression/prepare_eQTL.wdl", "eqtl_prepare_expression"): "standard",
    ("workflows/expression/rnaseqc2_aggregate_batched.wdl", "validate_rnaseqc_manifests"): "rnaseqc",
    ("workflows/expression/rnaseqc2_aggregate_batched.wdl", "aggregate_rnaseqc_batch"): "rnaseqc",
    ("workflows/expression/rnaseqc2_aggregate_batched.wdl", "merge_rnaseqc_batches"): "rnaseqc",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "BuildMethylationCohortSamples"): "standard",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "MergeMethylationChromosome"): "standard",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "BuildMethylationCorrelationCovariates"): "standard",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "AnalyzeMethylationCpGCorrelation"): "standard",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "AggregateMethylationChromosomes"): "standard",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "FinalizeMethylationConnectivity"): "standard",
    ("workflows/methylation/AggregateMethylationCohortArrays.wdl", "AnnotateMethylationSites"): "standard",
    ("workflows/methylation/ProcessMethylationSample.wdl", "FilterMethylationSample"): "methylation_rust",
    ("workflows/methylation/annotation.wdl", "AnnotateMethylationSites"): "standard",
    ("workflows/methylation/cohort_aggregation.wdl", "PrepareMethylationCohortManifest"): "standard",
    ("workflows/methylation/cohort_aggregation.wdl", "BuildMethylationCohortSamples"): "standard",
    ("workflows/methylation/cohort_aggregation.wdl", "MergeMethylationChromosome"): "methylation_rust",
    ("workflows/methylation/cohort_aggregation.wdl", "AggregateMethylationChromosomes"): "standard",
    ("workflows/methylation/connectivity.wdl", "BuildMethylationCorrelationCovariates"): "standard",
    ("workflows/methylation/connectivity.wdl", "AnalyzeMethylationCpGCorrelation"): "standard",
    ("workflows/methylation/connectivity.wdl", "FinalizeMethylationConnectivity"): "standard",
    ("workflows/methylation/merge_methylation.wdl", "ShardMethylationManifest"): "standard",
    ("workflows/methylation/merge_methylation.wdl", "FilterMethylationShard"): "standard",
    ("workflows/methylation/prepare_mQTL.wdl", "PrepareMethylationData"): "standard",
    ("workflows/proteomics/normalize_pQTL.wdl", "NormalizeProteomics"): "standard",
    ("workflows/proteomics/prepare_pQTL.wdl", "PrepareProteomicData"): "standard",
    ("workflows/splicing/prepare_sQTL.wdl", "PrepareSpliceData"): "standard",
}


def nested_calls(nodes):
    for node in nodes:
        if isinstance(node, WDL.Tree.Call):
            yield node
        elif isinstance(node, (WDL.Tree.Scatter, WDL.Tree.Conditional)):
            yield from nested_calls(node.body)


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

    def workflow_image_inputs(self, workflow):
        result = {}
        for declaration in workflow.inputs:
            value = literal_value(declaration)
            if value and value.startswith("ghcr.io/aou-multiomics-analysis/"):
                self.assertRegex(value, DIGEST_REFERENCE)
                result[declaration.name] = self.image_for_reference(value)
        return result

    def test_every_repository_runtime_uses_a_declared_string_input(self):
        actual = set()
        for path, document in self.documents.items():
            for task in document.tasks:
                if not task.runtime or "docker" not in task.runtime:
                    continue
                key = (path, task.name)
                actual.add(key)
                self.assertIn(key, TASK_IMAGES)
                image_input = str(task.runtime["docker"])
                declarations = {declaration.name: declaration for declaration in task.inputs}
                self.assertIn(image_input, declarations, key)
                self.assertEqual(str(declarations[image_input].type), "String", key)
                default = literal_value(declarations[image_input])
                if default is not None:
                    self.assertRegex(default, DIGEST_REFERENCE, key)
                    self.assertEqual(self.image_for_reference(default), TASK_IMAGES[key], key)
        self.assertEqual(actual, set(TASK_IMAGES))

    def test_calls_forward_the_matching_workflow_image(self):
        for path, document in self.documents.items():
            if document.workflow is None:
                continue
            # The cell workflows have four inputs for one repository. Their
            # stage-specific call selection is tested in test_stage_image_inputs.
            if path.startswith("workflows/cell_type_specific_expression/"):
                continue
            caller_images = self.workflow_image_inputs(document.workflow)
            environment = WDL.Env.Bindings()
            for name, image in caller_images.items():
                environment = environment.bind(name, WDL.Value.String("sentinel:" + name + ":" + image))
            for call in nested_calls(document.workflow.body):
                if isinstance(call.callee, WDL.Tree.Task):
                    callee_path = Path(call.callee.pos.abspath).relative_to(ROOT).as_posix()
                    key = (callee_path, call.callee.name)
                    if key not in TASK_IMAGES:
                        continue
                    wanted_image = TASK_IMAGES[key]
                    image_input = str(call.callee.runtime["docker"])
                    required = {image_input: wanted_image}
                else:
                    required = self.workflow_image_inputs(call.callee)
                for image_input, wanted_image in required.items():
                    matches = [name for name, image in caller_images.items() if image == wanted_image]
                    self.assertEqual(len(matches), 1, (path, call.name, image_input, wanted_image))
                    self.assertIn(image_input, call.inputs, (path, call.name))
                    actual = call.inputs[image_input].eval(environment, WDL.StdLib.Base("1.0")).value
                    self.assertEqual(actual, "sentinel:" + matches[0] + ":" + wanted_image,
                                     (path, call.name, image_input))

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
