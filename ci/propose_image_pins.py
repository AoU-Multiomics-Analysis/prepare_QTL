#!/usr/bin/env python3
"""Create candidate pin patches. No registry writes, Git commits, or WDL writes."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

import WDL
import yaml
from plan_image_updates import plan_changes

ROOT = Path(__file__).resolve().parents[1]


def validate_image(value, repository):
    if not isinstance(value, str) or not re.fullmatch(
            re.escape(repository) + r'@sha256:[0-9a-f]{64}', value):
        raise ValueError('Expected an immutable digest in repository ' + repository)


def literal_span(text, input_name):
    doc = WDL.parse_document(text)
    if doc.wdl_version != '1.0' or doc.workflow is None:
        raise ValueError('Pin targets must be WDL 1.0 workflows')
    declarations = [d for d in doc.workflow.inputs if d.name == input_name]
    if len(declarations) != 1:
        raise ValueError('Expected exactly one workflow input: ' + input_name)
    decl = declarations[0]
    if str(decl.type) != 'String' or decl.expr is None or not isinstance(decl.expr.literal, WDL.Value.String):
        raise ValueError('Pin target must have a literal String default: ' + input_name)
    pos = decl.expr.pos
    lines = text.splitlines(keepends=True)
    start = sum(map(len, lines[:pos.line - 1])) + pos.column - 1
    end = sum(map(len, lines[:pos.end_line - 1])) + pos.end_column - 1
    return start, end, decl.expr.literal.value


def propose(config, targets, plan, files, candidates):
    """Return changed file text only; caller-owned files and mappings are untouched."""
    if targets.get('version') != 1 or plan['unmapped']:
        raise ValueError('Invalid target version or unmapped source changes')
    unknown = set(plan['stages']) - set(targets['stages'])
    if unknown:
        raise ValueError('Stage pin targets are not configured: ' + ', '.join(sorted(unknown)))
    for image, value in candidates.items():
        if image not in config['images']:
            raise ValueError('Unknown candidate image: ' + image)
        validate_image(value, config['images'][image]['repository'])
    seen, edits = set(), {}
    for stage, locations in targets['stages'].items():
        if stage not in config['stages'] or not locations:
            raise ValueError('Unknown stage or empty pin target list: ' + stage)
        image = config['stages'][stage]['image']
        repository = config['images'][image]['repository']
        previous = set()
        for location in locations:
            path, name = location['path'], location['input']
            parts = PurePosixPath(path).parts
            if PurePosixPath(path).is_absolute() or '..' in parts or not path.endswith('.wdl'):
                raise ValueError('Invalid WDL target path: ' + path)
            if (path, name) in seen:
                raise ValueError('Duplicate pin target: ' + path + ':' + name)
            seen.add((path, name))
            start, end, old = literal_span(files[path], name)
            validate_image(old, repository)
            previous.add(old)
            if stage in plan['stages']:
                if image not in candidates:
                    raise ValueError('Missing candidate digest for ' + image)
                if candidates[image] != old:
                    edits.setdefault(path, []).append((start, end, json.dumps(candidates[image])))
        if len(previous) != 1:
            raise ValueError('Entry point defaults disagree for stage ' + stage)
    result = {}
    for path, spans in edits.items():
        updated = files[path]
        for start, end, replacement in sorted(spans, reverse=True):
            updated = updated[:start] + replacement + updated[end:]
        WDL.parse_document(updated)
        result[path] = updated
    return result


def git(repo, *args):
    return subprocess.check_output(['git', '-C', str(repo), *args], text=True).strip()


def checked_revisions(repo, base, head, expected_base, expected_head):
    resolved = tuple(git(repo, 'rev-parse', '--verify', '--end-of-options', ref + '^{commit}')
                     for ref in (base, head))
    if resolved != (expected_base, expected_head):
        raise ValueError('Base or head changed; discard the stale release proposal')
    subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor', *resolved], check=True)
    return resolved


def committed_text(repo, revision, path):
    return subprocess.check_output(['git', '-C', str(repo), 'show', revision + ':' + path], text=True)


def make_patch(files, updated):
    lines = []
    for path in sorted(updated):
        for line in difflib.unified_diff(files[path].splitlines(keepends=True),
                                        updated[path].splitlines(keepends=True),
                                        fromfile='a/' + path, tofile='b/' + path):
            lines.append(line if line.endswith('\n') else line + '\n\\ No newline at end of file\n')
    return ''.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT)
    parser.add_argument('--base', required=True)
    parser.add_argument('--head', default='HEAD')
    parser.add_argument('--expected-base', required=True)
    parser.add_argument('--expected-head', required=True)
    parser.add_argument('--candidate-image', action='append', default=[], metavar='IMAGE_ID=DIGEST_REFERENCE')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        base, head = checked_revisions(args.repo, args.base, args.head, args.expected_base, args.expected_head)
        candidates = {}
        for entry in args.candidate_image:
            image, separator, value = entry.partition('=')
            if not separator or not image or image in candidates:
                raise ValueError('Invalid or duplicate --candidate-image')
            candidates[image] = value
        config = yaml.safe_load(committed_text(args.repo, head, 'ci/image-stages.yml'))
        targets = yaml.safe_load(committed_text(args.repo, head, 'ci/release-pins.yml'))
        tracked = git(args.repo, 'ls-tree', '-r', '--name-only', '-z', head).split('\0')
        if plan_changes(config, [p for p in tracked if p])['unmapped']:
            raise ValueError('Committed source coverage is incomplete')
        changed = git(args.repo, 'diff', '--no-renames', '--name-only', '-z', base + '...' + head, '--').split('\0')
        plan = plan_changes(config, [p for p in changed if p])
        paths = {target['path'] for locations in targets['stages'].values() for target in locations}
        files = {path: committed_text(args.repo, head, path) for path in paths}
        updated = propose(config, targets, plan, files, candidates)
        patch = make_patch(files, updated)
        record = {'version': 1, 'status': 'candidate_not_validated', 'base': base, 'head': head,
                  'plan': plan, 'candidate_images': candidates, 'changed_files': sorted(updated),
                  'patch_sha256': hashlib.sha256(patch.encode()).hexdigest()}
        # Recheck local refs before emitting artifacts. Remote PR freshness is a later gate.
        checked_revisions(args.repo, args.base, args.head, base, head)
        args.output_dir.mkdir(parents=True, exist_ok=False)
        (args.output_dir / 'pins.patch').write_text(patch)
        (args.output_dir / 'release-candidate.json').write_text(json.dumps(record, indent=2) + '\n')
        print('Candidate files:', len(updated), '- no WDLs or GitHub state changed')
        return 0
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError, WDL.Error.SyntaxError) as error:
        print('Proposal rejected:', error, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
