"""Reuse the preceding smoke model and verify export/QTL outputs without a fit."""
import gzip
import json
from pathlib import Path
import subprocess

prefix = 'PrepareCellTypeEqtlWorkflow.'
baseline = json.loads(Path('ci-runs/precomputed/outputs.json').read_text())['outputs']
inputs = json.loads(Path('tests/cell_type_specific_expression/fixtures/precomputed-e2e.inputs.json').read_text())
inputs[prefix + 'precomputed_tca_model'] = baseline[prefix + 'tca_model']
# Both files exist, but are invalid as proportions/covariates. They must be ignored.
inputs[prefix + 'precomputed_proportions'] = inputs[prefix + 'lm22']
inputs[prefix + 'deconvolution_covariates'] = inputs[prefix + 'lm22']
inputs[prefix + 'gene_type'] = ['ignored_for_model_restart']
inputs[prefix + 'OutputPrefix'] = 'synthetic.restart'
input_path = Path('ci-runs/model-restart.inputs.json')
input_path.write_text(json.dumps(inputs, indent=2))
subprocess.run([
    'miniwdl', 'run', 'workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl',
    '--input', str(input_path), '--dir', 'ci-runs/model-restart/.', '--verbose', '--no-color'
], check=True)
result = json.loads(Path('ci-runs/model-restart/outputs.json').read_text())['outputs']
for name in ('estimated_proportions', 'tca_model_unfiltered', 'fit_tca_log',
             'proportions_lm22', 'proportions_combined', 'gene_type_filter_log'):
    assert result[prefix + name] is None, f'{name} should be skipped'
for name in ('cell_type_beds', 'filtered_cell_type_beds'):
    expected = baseline[prefix + name]
    actual = result[prefix + name]
    assert len(actual) == len(expected) > 0
    for left, right in zip(expected, actual):
        with gzip.open(left, 'rt') as a, gzip.open(right, 'rt') as b:
            assert a.read() == b.read(), f'Restarted {name} differs from the original export'
for name in ('int_beds', 'scaled_beds', 'int_phenotype_pcs', 'scaled_phenotype_pcs',
             'int_merged_covariates', 'scaled_merged_covariates'):
    assert len(result[prefix + name]) == len(result[prefix + 'cell_type_beds'])
    assert all(Path(path).is_file() for path in result[prefix + name])
assert Path(result[prefix + 'cell_type_qtl_manifest']).is_file()
parameters = json.loads(Path(result[prefix + 'effective_parameters_file']).read_text())
assert parameters['proportion_mode'] == 'precomputed_model'
assert 'tca_max_iters' not in parameters, 'Do not report new fitting settings for an old model'
print('Model restart: skipped fit; identical exported BEDs; QTL outputs and manifest present.')
