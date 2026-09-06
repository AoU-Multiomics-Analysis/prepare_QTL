#!/usr/bin/env python3
"""Run registered image checks against a patched, isolated source checkout."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import yaml
from plan_image_updates import matches
from propose_image_pins import literal_span, validate_image
from wdl_stage_routing import validate_routing

SUPPORTED_STAGES = {'cell_estimation', 'cell_fit', 'cell_export', 'cell_downstream',
                    'expression', 'common', 'proteomics', 'splicing', 'methylation',
                    'methylation_rust', 'rnaseqc'}

def selected_stages(record, config):
    stages = set(record['plan']['stages'])
    for path in record['plan']['wdl_checks']:
        for group in config['workflow_groups']:
            if matches(path, group['paths']):
                stages.update(group['stages'])
    return stages


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--suite', choices=('compact', 'full'),
                        help='Override cell-type coverage; releases default to compact, all-stages to full')
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--record', type=Path)
    selection.add_argument('--all-stages', action='store_true', help='Test all current defaults without publishing')
    args = parser.parse_args()
    trusted = Path(__file__).resolve().parents[1]
    source = args.source.resolve()
    routing_errors = validate_routing(source, trusted)
    if routing_errors:
        raise ValueError('Invalid stage routing:\n' + '\n'.join(routing_errors))
    config = yaml.safe_load((trusted / 'ci/image-stages.yml').read_text())
    targets = yaml.safe_load((trusted / 'ci/release-pins.yml').read_text())
    stages = set(config['stages']) if args.all_stages else selected_stages(json.loads(args.record.read_text()), config)
    if stages - SUPPORTED_STAGES:
        raise ValueError('No runtime gate for stages: ' + ', '.join(sorted(stages - SUPPORTED_STAGES)))
    for test in ('test_repo_image_routing.py', 'test_stage_image_inputs.py'):
        subprocess.run([sys.executable, str(trusted / 'tests' / test)], check=True,
                       env={**os.environ, 'RELEASE_SOURCE_ROOT': str(source)})
    # Policy tests are from trusted main. Candidate CI/tests cannot change in a release.
    for wdl in sorted((source / 'workflows').rglob('*.wdl')):
        subprocess.run(['miniwdl', 'check', str(wdl)], check=True)
    subprocess.run([sys.executable, str(trusted / 'scripts/check_wdl_file_scope.py'), str(source / 'workflows')], check=True)
    images = {}
    for stage in stages:
        values = set()
        for target in targets['stages'][stage]:
            value = literal_span((source / target['path']).read_text(), target['input'],
                                 target.get('scope', 'workflow'), target.get('task'))[2]
            validate_image(value, config['images'][config['stages'][stage]['image']]['repository'])
            values.add(value)
        if len(values) != 1:
            raise ValueError('Stage defaults disagree: ' + stage)
        images[stage] = values.pop()
    for image in set(images.values()):
        subprocess.run(['docker', 'pull', image], check=True)

    def run_r(image, script, *arguments, environment=None):
        command = ['docker', 'run', '--rm', '--user', f'{os.getuid()}:{os.getgid()}',
                   '--volume', f'{source}:{source}', '--workdir', str(source)]
        for key, value in (environment or {}).items():
            command += ['--env', key + '=' + value]
        subprocess.run(command + [image, 'Rscript', script, *arguments], check=True)

    if any(s.startswith('cell_') for s in stages) or stages & {'expression', 'common'}:
        # Runner reads each actual WDL digest default, not one image override.
        subprocess.run([sys.executable, str(trusted / 'tests/cell_type_specific_expression/smoke/run_pinned_images.py'),
                        '--suite', args.suite or ('full' if args.all_stages else 'compact')],
                       cwd=source, check=True)
    for image in sorted({images[s] for s in stages & {'expression', 'common'}}):
        for test in ('test_prepare_expression_log2_cpm.R', 'test_prepare_expression_sample_list.R'):
            run_r(image, 'tests/' + test,
                  environment={'PREPARE_EXPRESSION_SCRIPT': '/tmp/PrepareExpression.R'})
    if 'methylation' in stages:
        run_r(images['methylation'], 'tests/test_prepare_methylation.R',
              environment={'PREPARE_METHYLATION_SCRIPT': '/tmp/PrepareMethylation.R'})
    for stage in stages & {'proteomics', 'splicing'}:
        run_r(images[stage], 'tests/release/test_bundled_r.R', stage)
    if 'rnaseqc' in stages:
        subprocess.run([sys.executable, str(trusted / 'tests/rnaseqc2_aggregation/smoke_container.py'),
                        '--image', images['rnaseqc'], '--source-root', str(source)], cwd=source, check=True)
    if 'methylation_rust' in stages:
        subprocess.run([sys.executable, str(trusted / 'tests/release/test_rust_image.py'),
                        '--image', images['methylation_rust']], cwd=source, check=True)
    print('All selected image checks passed:', ', '.join(sorted(stages)))


if __name__ == '__main__':
    main()
