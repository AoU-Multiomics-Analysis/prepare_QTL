#!/usr/bin/env python3
"""Read-only image release impact report. Never build, publish, or edit WDLs."""
import argparse
from fnmatch import fnmatchcase
from pathlib import Path
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]


def matches(path, patterns):
    # Registry patterns use fnmatch semantics: * can span path separators.
    return any(fnmatchcase(path, pattern) for pattern in patterns)


def plan_changes(config, changed_paths):
    result = {key: set() for key in ('builds', 'stages', 'wdl_checks', 'test_changes', 'unmapped')}
    for path in sorted(set(changed_paths)):
        if path.endswith('.wdl'):
            result['wdl_checks'].add(path)
            if not any(matches(path, group['paths']) for group in config['workflow_groups']):
                result['unmapped'].add(path)
            continue
        if path.startswith('tests/'):
            result['test_changes'].add(path)
            continue
        classified = matches(path, config['ignored_sources'])
        for image, info in config['images'].items():
            if matches(path, info['build_paths']):
                result['builds'].add(image)
            if matches(path, info['environment_paths']):
                classified = True
                result['stages'].update(stage for stage, item in config['stages'].items()
                                        if item['image'] == image)
        for stage, info in config['stages'].items():
            if matches(path, info['sources']):
                classified = True
                result['stages'].add(stage)
        for rule in config['shared']:
            if matches(path, rule['sources']):
                classified = True
                result['stages'].update(rule['stages'])
        if path.startswith(('scripts/', 'rust/', 'envs/')) and not classified:
            result['unmapped'].add(path)
    return {key: sorted(value) for key, value in result.items()}


def validate_registry(config, root):
    errors = []
    if config.get('version') != 1:
        errors.append('Registry version must be 1')
    for name, stage in config['stages'].items():
        if stage['image'] not in config['images']:
            errors.append(f'Unknown image in stage {name}: {stage["image"]}')
    for group in config['shared'] + config['workflow_groups']:
        for stage in group['stages']:
            if stage not in config['stages']:
                errors.append(f'Unknown stage: {stage}')
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0')
    errors.extend('Unmapped source/workflow: ' + path for path in
                  plan_changes(config, tracked)['unmapped'])
    return errors


def report(plan):
    lines = ['## Image release plan — report only', '',
             'No images are built or published. No WDL pins are changed.', '',
             'Build dependencies and proposed pin updates are separate:', '']
    for key, title in [('builds', 'Image builds that existing source dependencies would require'),
                       ('stages', 'Release groups that may need new pins'),
                       ('wdl_checks', 'WDL changes requiring validation'),
                       ('test_changes', 'Changed tests'), ('unmapped', 'Unmapped paths — review required')]:
        lines.extend([f'### {title}', ''])
        lines.extend('- `' + item.replace('`', '\\`') + '`' for item in plan[key])
        if not plan[key]:
            lines.append('None.')
        lines.append('')
    lines.append('This is a dependency proposal, not proof of cache reuse or mixed-image compatibility.')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'ci/image-stages.yml')
    parser.add_argument('--base')
    parser.add_argument('--head', default='HEAD')
    parser.add_argument('--validate', action='store_true')
    parser.add_argument('--changed', action='append', default=[], help='Simulate one changed path; repeatable')
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    errors = validate_registry(config, ROOT)
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        return 1
    changed = args.changed
    if args.base:
        base = subprocess.check_output(['git', 'rev-parse', '--verify', '--end-of-options', args.base + '^{commit}'], cwd=ROOT).decode().strip()
        head = subprocess.check_output(['git', 'rev-parse', '--verify', '--end-of-options', args.head + '^{commit}'], cwd=ROOT).decode().strip()
        changed += subprocess.check_output(['git', 'diff', '--no-renames', '--name-only', '-z', base + '...' + head, '--'], cwd=ROOT).decode().split('\0')
    if not (args.base or args.changed or args.validate):
        parser.error('Provide --base, --changed, or --validate')
    plan = plan_changes(config, [path for path in changed if path])
    print(report(plan))
    return int(bool(plan['unmapped']))


if __name__ == '__main__':
    sys.exit(main())
