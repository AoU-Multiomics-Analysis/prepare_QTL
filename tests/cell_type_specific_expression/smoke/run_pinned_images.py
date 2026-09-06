"""GitHub-only smoke: pull WDL digest defaults, never build or replace images."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import WDL

PREFIX = 'PrepareCellTypeEqtlWorkflow.'
WORKFLOW = 'workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl'
DECONVOLUTION = 'workflows/cell_type_specific_expression/deconvolution.wdl'
DECONV_PREFIX = 'CellTypeDeconvolution.'
RESTART_RUNNER = 'tests/cell_type_specific_expression/smoke/run_model_restart.py'


def pinned_inputs(fixture, pins):
    return {**fixture, **{PREFIX + key: value for key, value in pins.items()}}


def compact_inputs(fixture, mode):
    inputs = dict(fixture)
    for key in inputs:
        if key.endswith('_cpu'):
            inputs[key] = 1
        elif key.endswith('_memory'):
            inputs[key] = 4 if isinstance(inputs[key], int) else '4 GB'
    inputs[PREFIX + 'gene_summary_memory'] = '4 GB'
    if mode == 'precomputed':
        # Fixed fixture cohort means retain B, CD4, and monocyte/myeloid groups.
        inputs[PREFIX + 'group_mean_threshold'] = .12
    return inputs


def scoped_inputs(inputs, workflow, prefix):
    allowed = {d.name for d in WDL.load(workflow).workflow.inputs}
    result = {}
    for key, value in inputs.items():
        name = key.removeprefix(PREFIX)
        if name == 'deconvolution_covariates':
            name = 'covariates'
        if name in allowed:
            result[prefix + name] = value
    return result


def read_outputs(path):
    # MiniWDL's outputs.json is flat; only its CLI stdout has an outputs wrapper.
    return json.loads(Path(path).read_text())


def restart_command(baseline_outputs, baseline_inputs, restart_inputs, output_directory,
                    deconvolution_only=False):
    command = [
        sys.executable,
        RESTART_RUNNER,
        '--baseline-outputs', str(baseline_outputs),
        '--baseline-inputs', str(baseline_inputs),
        '--restart-inputs', str(restart_inputs),
        '--output-directory', str(output_directory),
    ]
    return command + (['--deconvolution-only'] if deconvolution_only else [])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('compact', 'full'), default='full')
    args = parser.parse_args()
    print(f'stage=smoke_plan suite={args.suite} status=start', flush=True)
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
        fixture_directory = Path('tests/cell_type_specific_expression/fixtures')
        if args.suite == 'compact':
            inputs = compact_inputs(inputs, mode)
            if mode == 'precomputed':
                fixture_directory = Path('ci-runs/compact-fixture')
                shutil.copytree('tests/cell_type_specific_expression/fixtures', fixture_directory,
                                dirs_exist_ok=True)
                (fixture_directory / 'expected_groups.txt').write_text('B cells\nCD4 T cells\nMonocyte/myeloid\n')
        input_path = Path(f'ci-runs/pinned-{mode}.inputs.json')
        input_path.write_text(json.dumps(inputs, indent=2))
        if mode == 'precomputed':
            run_r(pins['downstream_docker_image'],
                  'tests/cell_type_specific_expression/generate_reference_fixture.R',
                  'tests/cell_type_specific_expression/fixtures/synthetic_expression.bed',
                  'ci-runs/reference-counts.tsv.gz')
        run_dir = f'ci-runs/pinned-{mode}'
        workflow = WORKFLOW
        run_inputs = input_path
        if args.suite == 'compact' and mode == 'hspe':
            workflow = DECONVOLUTION
            run_inputs = Path('ci-runs/pinned-hspe.deconvolution.inputs.json')
            run_inputs.write_text(json.dumps(scoped_inputs(inputs, workflow, DECONV_PREFIX), indent=2))
        print(f'stage=smoke_workflow mode={mode} workflow={workflow} suite={args.suite}', flush=True)
        subprocess.run(['miniwdl', 'run', workflow, '--input', str(run_inputs),
                        '--dir', run_dir + '/.', '--verbose', '--no-color'], check=True)
        output_path = run_dir + '/outputs.json'
        outputs = read_outputs(output_path)
        expected = {key.removesuffix('_docker_image'): value for key, value in pins.items()}
        assertions = [('assert_deconvolution_outputs.R', pins['downstream_docker_image']),
                      ('assert_qtl_outputs.R', pins['qtl_docker_image'])]
        if workflow == DECONVOLUTION:
            outputs = {PREFIX + key.removeprefix(DECONV_PREFIX): value for key, value in outputs.items()}
            output_path = run_dir + '/assertion-outputs.json'
            Path(output_path).write_text(json.dumps(outputs, indent=2))
            expected.pop('qtl')
            assertions = assertions[:1]
        if outputs[PREFIX + 'stage_images'] != expected:
            raise AssertionError('Workflow did not retain the selected stage images')
        for script, image in assertions:
            run_r(image, 'tests/cell_type_specific_expression/smoke/' + script,
                  output_path, str(input_path), str(fixture_directory))
        if mode == 'precomputed':
            subprocess.run(restart_command(
                Path(output_path),
                input_path,
                Path('ci-runs/pinned-model-restart.inputs.json'),
                Path('ci-runs/pinned-model-restart'),
                deconvolution_only=args.suite == 'compact',
            ), check=True)
        print(f'stage=pinned_smoke mode={mode} status=passed', flush=True)


if __name__ == '__main__':
    main()
