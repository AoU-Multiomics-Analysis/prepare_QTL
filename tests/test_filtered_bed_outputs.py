"""Filtered BED discovery and scatter retain actual cloud File identities."""
from pathlib import Path
import tempfile
import unittest

import WDL

ROOT = Path(__file__).resolve().parents[1]


class OutputStdLib(WDL.StdLib.TaskOutputs):
    def __init__(self, directory):
        super().__init__('1.0')
        self.directory = Path(directory)
        self.glob = WDL.StdLib.StaticFunction('glob', [WDL.Type.String()],
                                             WDL.Type.Array(WDL.Type.File()), self.discover)

    def _devirtualize_filename(self, filename):
        return str(self.directory / filename)

    def discover(self, pattern):
        return WDL.Value.Array(WDL.Type.File(), [
            WDL.Value.File('gs://bucket/call-FilterCellTypeBeds/glob-actual/' + path.name)
            for path in sorted(self.directory.glob(pattern.value))])


class FilteredBedOutputsTest(unittest.TestCase):
    def test_discover_files_instead_of_reconstructing_paths_from_text(self):
        doc = WDL.load(str(ROOT / 'workflows/cell_type_specific_expression/tasks/reference_filter.wdl'))
        task = next(t for t in doc.tasks if t.name == 'FilterCellTypeBeds')
        output = next(d for d in task.outputs if d.name == 'filtered_beds')
        with tempfile.TemporaryDirectory() as directory:
            beds = Path(directory) / 'outputs/beds'
            beds.mkdir(parents=True)
            for name in ('eosinophils', 'cd4_t_cells'):
                (beds / f'{name}.filtered.bed.gz').write_bytes(b'fixture')
            (Path(directory) / 'outputs/filtered_beds.txt').write_text('outputs/beds/not-uploaded.bed.gz\n')
            value = output.expr.eval(WDL.Env.Bindings(), OutputStdLib(directory)).coerce(output.type)
            self.assertIsInstance(value.type.item_type, WDL.Type.File)
            self.assertEqual(value.json, [
                'gs://bucket/call-FilterCellTypeBeds/glob-actual/cd4_t_cells.filtered.bed.gz',
                'gs://bucket/call-FilterCellTypeBeds/glob-actual/eosinophils.filtered.bed.gz'])

    def test_match_files_to_labels_independent_of_glob_order(self):
        doc = WDL.load(str(ROOT / 'workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl'))
        outer = next(n for n in doc.workflow.body if isinstance(n, WDL.Tree.Scatter))
        inner = next((n for n in outer.body if isinstance(n, WDL.Tree.Scatter)), None)
        self.assertIsNotNone(inner, 'BEDs must be matched by identity, not array position')
        condition = next(n for n in inner.body if isinstance(n, WDL.Tree.Conditional))
        call = next(n for n in outer.body if isinstance(n, WDL.Tree.Call))
        slugs = ['eosinophils', 'cd4_t_cells']  # inventory differs from glob order
        files = ['gs://bucket/glob-actual/cd4_t_cells.filtered.bed.gz',
                 'gs://bucket/glob-actual/eosinophils.filtered.bed.gz']
        stdlib = WDL.StdLib.Base('1.0')
        for index, slug in enumerate(slugs):
            env = WDL.Env.Bindings().bind('index', WDL.Value.Int(index)).bind(
                'PrepareScatterInputs.cell_type_slugs', WDL.Value.Array(
                    WDL.Type.String(), [WDL.Value.String(s) for s in slugs]))
            matches = []
            for path in files:
                local = env.bind('filtered_bed', WDL.Value.File(path))
                matches.append(WDL.Value.File(path) if condition.expr.eval(local, stdlib).value
                               else WDL.Value.Null())
            env = env.bind('matched_filtered_bed', WDL.Value.Array(WDL.Type.File(optional=True), matches))
            selected = call.inputs['CpmBed'].eval(env, stdlib)
            self.assertEqual(str(selected.type), 'File')
            self.assertEqual(selected.value, f'gs://bucket/glob-actual/{slug}.filtered.bed.gz')
            with tempfile.TemporaryDirectory() as directory:
                localized_bed = Path(directory) / f"donor's {slug}.filtered.bed.gz"
                localized_bed.write_bytes(b'localized fixture')
                localized = WDL.Value.rewrite_paths(selected, lambda value: str(localized_bed))
                self.assertEqual(localized.value, str(localized_bed))
                self.assertTrue(Path(localized.value).is_file())
        missing = WDL.Env.Bindings().bind('matched_filtered_bed', WDL.Value.Array(
            WDL.Type.File(optional=True), [WDL.Value.Null(), WDL.Value.Null()]))
        with self.assertRaises(WDL.Error.EvalError):
            call.inputs['CpmBed'].eval(missing, stdlib)


if __name__ == '__main__':
    unittest.main()
