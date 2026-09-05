"""GitHub-only smoke: pull WDL digest defaults, never build or replace images."""
import json
import os
from pathlib import Path
import re
import subprocess
import WDL

PREFIX = 'PrepareCellTypeEqtlWorkflow.'
WORKFLOW = 'workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl'


def pinned_inputs(fixture, pins):
    return {**fixture, **{PREFIX + key: value for key, value in pins.items()}}


def read_outputs(path):
    # MiniWDL's outputs.json is flat; only its CLI stdout has an outputs wrapper.
    return json.loads(Path(path).read_text())


def main():
    root = Path.cwd()
    doc = WDL.load(WORKFLOW)
    pins = {decl.name: decl.expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value
            for decl in doc.workflow.inputs if decl.name.endswith('_docker_image')}
    for image in set(pins.values()):
        if not re.fullmatch(r'ghcr\.io/[^\s@]+@sha256:[0-9a-f]{64}', image):
            raise ValueError('Smoke requires immutable image defaults: ' + image)
        subprocess.run(['docker', 'pull', image], check=True)

    def run_r(image, script, *args):
        subprocess.run(['docker', 'run', '--rm', '--user', f'{os.getuid()}:{os.getgid()}',
                        '--group-add', str(os.getgid()), '--volume', f'{root}:{root}',
                        '--workdir', str(root), image, 'Rscript', script, *args], check=True)

    Path('ci-runs').mkdir(exist_ok=True)
    for mode in ('hspe', 'precomputed'):
        print(f'stage=pinned_smoke mode={mode} status=start', flush=True)
        fixture = Path(f'tests/cell_type_specific_expression/fixtures/{mode}-e2e.inputs.json')
        inputs = pinned_inputs(json.loads(fixture.read_text()), pins)
        input_path = Path(f'ci-runs/pinned-{mode}.inputs.json')
        input_path.write_text(json.dumps(inputs, indent=2))
        if mode == 'precomputed':
            run_r(pins['downstream_docker_image'],
                  'tests/cell_type_specific_expression/generate_reference_fixture.R',
                  'tests/cell_type_specific_expression/fixtures/synthetic_expression.bed',
                  'ci-runs/reference-counts.tsv.gz')
        run_dir = f'ci-runs/pinned-{mode}'
        subprocess.run(['miniwdl', 'run', WORKFLOW, '--input', str(input_path),
                        '--dir', run_dir + '/.', '--verbose', '--no-color'], check=True)
        output_path = run_dir + '/outputs.json'
        outputs = read_outputs(output_path)
        expected = {key.removesuffix('_docker_image'): value for key, value in pins.items()}
        if outputs[PREFIX + 'stage_images'] != expected:
            raise AssertionError('Workflow did not retain the selected stage images')
        for script, image in [('assert_deconvolution_outputs.R', pins['downstream_docker_image']),
                              ('assert_qtl_outputs.R', pins['qtl_docker_image'])]:
            run_r(image, 'tests/cell_type_specific_expression/smoke/' + script,
                  output_path, str(input_path), 'tests/cell_type_specific_expression/fixtures')
        print(f'stage=pinned_smoke mode={mode} status=passed', flush=True)


if __name__ == '__main__':
    main()
