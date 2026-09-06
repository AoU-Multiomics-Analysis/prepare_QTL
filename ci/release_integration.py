#!/usr/bin/env python3
"""Trusted release controller. Run this from main, never from the candidate tree."""
import argparse
import base64
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.error
import urllib.request

import yaml
from plan_image_updates import matches, plan_changes
from propose_image_pins import committed_text, git, make_patch, propose, validate_image


def validate_pr(pr, repository, previous=None):
    if (pr['state'] != 'open' or pr['base']['ref'] != 'main'
            or pr['base']['repo']['full_name'].lower() != repository.lower()
            or pr['head']['repo']['full_name'].lower() != repository.lower()
            or pr['head']['ref'] == 'main'):
        raise ValueError('Release requires an open same-repository PR targeting main')
    for side in ('base', 'head'):
        if not re.fullmatch(r'[0-9a-f]{40}', pr[side]['sha']):
            raise ValueError('Invalid PR commit')
        if previous and (pr[side]['sha'] != previous[side]['sha'] or pr[side]['ref'] != previous[side]['ref']):
            raise ValueError('Stale release: PR head or base changed')


def check_policy_changes(paths):
    blocked = [p for p in paths if p.startswith(('ci/', '.github/', 'tests/'))]
    if blocked:
        raise ValueError('Merge release policy/test changes separately before releasing source: ' + ', '.join(blocked))


def source_fingerprint(repo, revision, patterns):
    entries = subprocess.check_output(['git', '-C', str(repo), 'ls-tree', '-r', '-z', revision])
    selected = []
    for entry in entries.split(b'\0'):
        if entry and matches(entry.split(b'\t', 1)[1].decode(), patterns):
            selected.append(entry)
    if not selected:
        raise ValueError('Image has no tracked build inputs')
    return hashlib.sha256(b'\0'.join(sorted(selected))).hexdigest()


def validate_build_result(spec, result):
    if result['fingerprint'] != spec['fingerprint']:
        raise ValueError('Image source fingerprint does not match this release')
    validate_image(result['reference'], spec['repository'])


def commit_payload(repository, branch, head, files):
    return {'branch': {'repositoryNameWithOwner': repository, 'branchName': branch},
            'expectedHeadOid': head, 'message': {'headline': 'Pin tested stage images'},
            'fileChanges': {'additions': [{'path': p, 'contents': base64.b64encode(text.encode()).decode()}
                                        for p, text in sorted(files.items())]}}


def api(repository, route):
    return json.loads(subprocess.check_output(['gh', 'api', f'repos/{repository}/{route}'], text=True))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def snapshot(repo, repository, number, output, expected_head=None):
    pr = api(repository, f'pulls/{number}')
    validate_pr(pr, repository)
    if expected_head and pr['head']['sha'] != expected_head:
        print('Release request superseded: skip queued build and tests')
        if os.environ.get('GITHUB_OUTPUT'):
            with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
                stream.write('skip=true\nhas_builds=false\n')
        return
    subprocess.run(['git', '-C', str(repo), 'fetch', 'origin', pr['base']['sha'], pr['head']['sha']], check=True)
    # The checked-out policy must be exactly the current PR base on main.
    if git(repo, 'rev-parse', 'HEAD') != pr['base']['sha']:
        raise ValueError('Trusted checkout is stale; dispatch again from current main')
    subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor', pr['base']['sha'], pr['head']['sha']], check=True)
    changed = git(repo, 'diff', '--no-renames', '--name-only', '-z', pr['base']['sha'], pr['head']['sha'], '--').split('\0')
    check_policy_changes(changed)
    config = yaml.safe_load((repo / 'ci/image-stages.yml').read_text())
    targets = yaml.safe_load((repo / 'ci/release-pins.yml').read_text())
    plan = plan_changes(config, [p for p in changed if p])
    if plan['unmapped'] or set(plan['stages']) - set(targets['stages']):
        raise ValueError('Unmapped source or missing stage targets')
    # Rebuild only images needed by changed stage consumers, not unused copies.
    wanted = sorted({config['stages'][s]['image'] for s in plan['stages']})
    builds = []
    for name in wanted:
        spec = config['images'][name]
        fingerprint = source_fingerprint(repo, pr['head']['sha'], spec['build_paths'])
        builds.append({'id': name, 'repository': spec['repository'], 'dockerfile': spec['dockerfile'],
                       'fingerprint': fingerprint, 'tag': spec['repository'] + ':release-src-' + fingerprint})
    result = {'version': 1, 'repository': repository, 'pr': pr, 'plan': plan, 'builds': builds}
    save(output, result)
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as stream:
            stream.write('head=' + pr['head']['sha'] + '\n')
            stream.write('base=' + pr['base']['sha'] + '\n')
            stream.write('matrix=' + json.dumps({'include': builds}, separators=(',', ':')) + '\n')
            stream.write('has_builds=' + str(bool(builds)).lower() + '\n')


def registry_image(reference):
    """Read public GHCR manifests and labels; absence is not a valid cached build."""
    name, tag = reference.removeprefix('ghcr.io/').rsplit(':', 1) if '@' not in reference else reference.removeprefix('ghcr.io/').split('@', 1)
    token = json.load(urllib.request.urlopen('https://ghcr.io/token?scope=repository:' + name + ':pull', timeout=30))['token']
    def read(path):
        request = urllib.request.Request('https://ghcr.io/v2/' + name + '/' + path,
            headers={'Authorization': 'Bearer ' + token, 'Accept': 'application/vnd.docker.distribution.manifest.v2+json, application/vnd.oci.image.manifest.v1+json'})
        with urllib.request.urlopen(request, timeout=30) as response:
            data = response.read()
        return json.loads(data), 'sha256:' + hashlib.sha256(data).hexdigest()
    manifest, digest = read('manifests/' + tag)
    config, _ = read('blobs/' + manifest['config']['digest'])
    return 'ghcr.io/' + name + '@' + digest, config.get('config', {}).get('Labels', {})


def build(source, spec, head, output):
    try:
        reference, labels = registry_image(spec['tag'])
        if labels.get('org.prepare-qtl.source-fingerprint') != spec['fingerprint']:
            raise ValueError('Existing source tag has unexpected provenance')
    except urllib.error.HTTPError as error:
        if error.code != 404:
            raise
        metadata = Path(output).with_suffix('.build-metadata.json')
        subprocess.run(['docker', 'buildx', 'build', '--platform', 'linux/amd64', '--provenance=false',
                        '--file', str(source / spec['dockerfile']), '--tag', spec['tag'], '--push',
                        '--label', 'org.prepare-qtl.source-fingerprint=' + spec['fingerprint'],
                        '--label', 'org.opencontainers.image.revision=' + head,
                        '--metadata-file', str(metadata), str(source)], check=True)
        digest = json.loads(metadata.read_text())['containerimage.digest']
        reference = spec['repository'] + '@' + digest
        resolved, labels = registry_image(reference)
        if resolved != reference or labels.get('org.prepare-qtl.source-fingerprint') != spec['fingerprint']:
            raise ValueError('Published image does not match the selected source')
    result = {'id': spec['id'], 'fingerprint': spec['fingerprint'], 'reference': reference}
    validate_build_result(spec, result)
    save(output, result)


def candidate(repo, record, build_directory):
    config = yaml.safe_load((repo / 'ci/image-stages.yml').read_text())
    targets = yaml.safe_load((repo / 'ci/release-pins.yml').read_text())
    candidates = {}
    for spec in record['builds']:
        result = json.loads((build_directory / (spec['id'] + '.json')).read_text())
        validate_build_result(spec, result)
        # Independently resolve immutable digest and source label, even in the updater.
        resolved, labels = registry_image(result['reference'])
        if resolved != result['reference'] or labels.get('org.prepare-qtl.source-fingerprint') != spec['fingerprint']:
            raise ValueError('Registry provenance check failed')
        candidates[spec['id']] = result['reference']
    paths = {target['path'] for locations in targets['stages'].values() for target in locations}
    files = {path: committed_text(repo, record['pr']['head']['sha'], path) for path in paths}
    updated = propose(config, targets, record['plan'], files, candidates)
    return files, updated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['snapshot', 'build', 'prepare', 'commit'])
    parser.add_argument('--repo', type=Path, default=Path.cwd())
    parser.add_argument('--repository', default=os.environ.get('GITHUB_REPOSITORY'))
    parser.add_argument('--pr', type=int)
    parser.add_argument('--expected-head', default='')
    parser.add_argument('--record', type=Path, default=Path('release.json'))
    parser.add_argument('--builds', type=Path, default=Path('build-results'))
    parser.add_argument('--source', type=Path)
    parser.add_argument('--image-id')
    parser.add_argument('--output', type=Path, default=Path('candidate'))
    args = parser.parse_args()
    if args.operation == 'snapshot':
        if not args.pr or args.pr <= 0:
            parser.error('--pr must be a positive PR number')
        snapshot(args.repo, args.repository, args.pr, args.record, args.expected_head)
        return
    record = json.loads(args.record.read_text())
    validate_pr(record['pr'], args.repository)
    if record['repository'] != args.repository:
        raise ValueError('Release repository mismatch')
    if args.operation == 'build':
        spec = next(item for item in record['builds'] if item['id'] == args.image_id)
        if git(args.source, 'rev-parse', 'HEAD') != record['pr']['head']['sha']:
            raise ValueError('Build checkout mismatch')
        build(args.source, spec, record['pr']['head']['sha'], args.output)
        return
    # Record is data: independently reproduce the plan from trusted base policy.
    fresh = args.record.with_suffix('.rechecked.json')
    snapshot(args.repo, args.repository, record['pr']['number'], fresh)
    checked = json.loads(fresh.read_text())
    validate_pr(checked['pr'], args.repository, record['pr'])
    if checked['plan'] != record['plan'] or checked['builds'] != record['builds']:
        raise ValueError('Release artifact disagrees with trusted policy')
    files, updated = candidate(args.repo, record, args.builds)
    if args.operation == 'prepare':
        args.output.mkdir(exist_ok=False)
        (args.output / 'pins.patch').write_text(make_patch(files, updated))
        save(args.output / 'images.json', {item['id']: json.loads((args.builds / (item['id'] + '.json')).read_text())['reference'] for item in record['builds']})
        return
    if os.environ.get('RELEASE_AUTO_COMMIT') != 'true':
        raise ValueError('Automatic commits are disabled')
    if not updated:
        print('Pins already current; no commit needed')
        return
    validate_pr(api(args.repository, f'pulls/{record["pr"]["number"]}'), args.repository, record['pr'])
    payload = commit_payload(args.repository, record['pr']['head']['ref'], record['pr']['head']['sha'], updated)
    query = 'mutation($input: CreateCommitOnBranchInput!) { createCommitOnBranch(input: $input) { commit { oid url } } }'
    response = json.loads(subprocess.check_output(['gh', 'api', 'graphql', '--input', '-'],
        input=json.dumps({'query': query, 'variables': {'input': payload}}), text=True))
    if response.get('errors'):
        raise ValueError('Pin commit rejected: ' + json.dumps(response['errors']))
    print(response['data']['createCommitOnBranch']['commit']['url'])


if __name__ == '__main__':
    main()
