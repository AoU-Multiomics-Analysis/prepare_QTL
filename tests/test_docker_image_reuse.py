import importlib.util
from pathlib import Path
import subprocess
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
IMAGE = 'ghcr.io/example/image@sha256:' + 'a' * 64


class ImageReuse(unittest.TestCase):
    def load_helper(self):
        path = ROOT / 'ci/docker_images.py'
        self.assertTrue(path.exists(), 'Shared pinned-image reuse helper is missing')
        spec = importlib.util.spec_from_file_location('image_reuse', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_second_runner_reuses_exact_image(self):
        helper = self.load_helper()
        present = set()
        commands = []

        def docker(command, **kwargs):
            commands.append(command)
            if command[:3] == ['docker', 'image', 'inspect']:
                return subprocess.CompletedProcess(command, 0 if command[-1] in present else 1)
            self.assertEqual(command[:2], ['docker', 'pull'])
            self.assertTrue(kwargs['check'])
            present.add(command[-1])
            return subprocess.CompletedProcess(command, 0)

        with mock.patch.object(helper.subprocess, 'run', side_effect=docker):
            helper.ensure_pinned_image(IMAGE)
            # No process-global cache: a separate runner can inspect Docker too.
            self.load_helper().ensure_pinned_image(IMAGE)
            helper.ensure_pinned_image(IMAGE.replace('a' * 64, 'b' * 64))
        self.assertEqual(sum(cmd[:2] == ['docker', 'pull'] for cmd in commands), 2)
        self.assertEqual(commands[0], ['docker', 'image', 'inspect', IMAGE])

    def test_tag_is_rejected_without_contacting_docker(self):
        helper = self.load_helper()
        with mock.patch.object(helper.subprocess, 'run') as execute:
            with self.assertRaisesRegex(ValueError, 'digest'):
                helper.ensure_pinned_image('ghcr.io/example/image:main')
            execute.assert_not_called()

    def test_pull_failure_stops_the_gate(self):
        helper = self.load_helper()
        with mock.patch.object(helper.subprocess, 'run', side_effect=[
            subprocess.CompletedProcess([], 1), subprocess.CalledProcessError(1, ['docker', 'pull'])
        ]):
            with self.assertRaises(subprocess.CalledProcessError):
                helper.ensure_pinned_image(IMAGE)


if __name__ == '__main__':
    unittest.main()
