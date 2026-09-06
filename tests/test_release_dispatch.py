"""Execute the actual metadata-only dispatcher with a fake GitHub transport."""
import json
from pathlib import Path
import subprocess
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
BOT = 'aou-prepare-qtl-release[bot]'


class DispatchTests(unittest.TestCase):
    def dispatch(self, *, action='synchronize', filename='workflows/example.wdl',
                 patch=None, author=BOT, parent='before', file_status='modified', additions=1):
        flow = yaml.load((ROOT / '.github/workflows/image-release-dispatch.yml').read_text(), Loader=yaml.BaseLoader)
        script = flow['jobs']['dispatch']['steps'][0]['with']['script']
        if patch is None:
            patch = ('@@ -1 +1 @@\n- String downstream_docker_image = "ghcr.io/org/image@sha256:' +
                     'a' * 64 + '"\n+ String downstream_docker_image = "ghcr.io/org/image@sha256:' + 'b' * 64 + '"')
        case = {'script': script, 'context': {'repo': {'owner': 'org', 'repo': 'repo'},
                'payload': {'action': action, 'before': 'before', 'after': 'head',
                            'sender': {'login': author},
                            'pull_request': {'number': 42, 'head': {'sha': 'head'}}}},
                'commit': {'sha': 'head', 'author': {'login': author},
                           'commit': {'message': 'Pin published stage images'},
                           'parents': [{'sha': parent}],
                           'files': [{'filename': filename, 'status': file_status, 'patch': patch,
                                      'additions': additions, 'deletions': 1}]}}
        harness = r'''
const fs = require('fs');
const x = JSON.parse(fs.readFileSync(0, 'utf8'));
const sent = [];
const github = {rest: {repos: {getCommit: async () => ({data: x.commit})},
  actions: {createWorkflowDispatch: async input => sent.push(input)}}};
const core = {info: () => {}};
const AsyncFunction = Object.getPrototypeOf(async function(){}).constructor;
new AsyncFunction('github', 'context', 'core', 'process', x.script)(
  github, x.context, core, {env: {COMMIT_PINS: 'true'}}
).then(() => process.stdout.write(JSON.stringify(sent))).catch(e => {throw e});
'''
        result = subprocess.run(['node', '-e', harness], input=json.dumps(case),
                                text=True, capture_output=True, check=True)
        return json.loads(result.stdout)

    def test_app_pin_only_synchronization_does_not_dispatch(self):
        self.assertEqual(self.dispatch(), [])

    def test_source_edit_still_dispatches(self):
        self.assertEqual(len(self.dispatch(filename='scripts/tool.R', patch='@@\n-a\n+b')), 1)

    def test_manual_pin_change_still_dispatches(self):
        self.assertEqual(len(self.dispatch(author='contributor')), 1)

    def test_explicit_label_can_retest_pins(self):
        self.assertEqual(len(self.dispatch(action='labeled')), 1)

    def test_multiple_pushed_commits_cannot_hide_source_changes(self):
        self.assertEqual(len(self.dispatch(parent='intermediate')), 1)

    def test_other_wdl_edits_still_dispatch(self):
        self.assertEqual(len(self.dispatch(patch='@@\n- Int cores = 1\n+ Int cores = 2')), 1)

    def test_incomplete_patch_and_renames_still_dispatch(self):
        self.assertEqual(len(self.dispatch(patch='')), 1)
        self.assertEqual(len(self.dispatch(file_status='renamed')), 1)

    def test_truncated_pin_patch_cannot_hide_other_changes(self):
        self.assertEqual(len(self.dispatch(additions=2)), 1)


if __name__ == '__main__':
    unittest.main()
