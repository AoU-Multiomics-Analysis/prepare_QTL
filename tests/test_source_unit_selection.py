"""Exercise the workflow's real source-selection shell against a local Git repo."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]


class SourceSelectionTests(unittest.TestCase):
    def test_source_runs_but_following_pin_commit_does_not(self):
        workflow = yaml.load((ROOT / '.github/workflows/source-unit-checks.yml').read_text(), Loader=yaml.BaseLoader)
        command = next(step['run'] for step in workflow['jobs']['unit']['steps']
                       if step.get('id') == 'changes')
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)

            def git(*args):
                return subprocess.check_output(['git', '-C', directory, *args], text=True).strip()

            git('init', '-q')
            git('config', 'user.name', 'Test')
            git('config', 'user.email', 'test@example.invalid')

            def commit(path, content):
                target = repo / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(content)
                git('add', path)
                git('commit', '-qm', 'fixture')
                return git('rev-parse', 'HEAD')

            base = commit('README', 'fixture')
            source = commit('scripts/cell_type_specific_expression/run_hspe.R', '# source')
            pins = commit('workflows/main.wdl', '# pins')
            for before, after, expected in ((base, source, 'true'), (source, pins, 'false'),
                                             (base, pins, 'true'), ('f' * 40, pins, 'true')):
                with self.subTest(before=before, after=after):
                    output = repo / 'job-output'
                    output.write_text('')
                    subprocess.run(['bash', '-c', command], cwd=repo, check=True,
                                   env={**os.environ, 'BASE': before, 'HEAD': after,
                                        'EVENT': 'pull_request', 'GITHUB_OUTPUT': str(output)})
                    self.assertEqual(output.read_text().strip(), f'run={expected}')


if __name__ == '__main__':
    unittest.main()
