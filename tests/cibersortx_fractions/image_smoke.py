"""Run the rendered wrapper in the official image on GitHub Actions.

Only the authenticated native calculation is replaced with a small fixture.
No image is built and no credentials or scientific data are required.
"""
import argparse
import importlib.util
import json
import os
import subprocess
import tempfile
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]
IMAGE = 'cibersortx/fractions@sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc'


def prepare(base, smode):
    document = WDL.load(str(ROOT / 'workflows/cell_type_specific_expression/cibersortx_fractions.wdl'))
    task = document.tasks[0]
    localized = base / 'localized'
    localized.mkdir()
    contents = {
        'mixture': 'GeneSymbol\ts1\ts2\ts3\nsig\t1\t2\t3\nZNF804A\t2\t3\t4\n',
        'signature': 'GeneSymbol\tB\tCD4_T\nsig\t1\t2\n',
        'refsample': 'GeneSymbol\tB\tB.1\tCD4_T\tCD4_T.1\nsig\t1\t2\t3\t4\n',
        'source_geps': 'GeneSymbol\tB\tCD4_T\nsig\t1\t2\n',
    }
    values = dict(username='fixture@example.org', token='fixture-token', smode=smode)
    for name, text in contents.items():
        (localized / (name + '.txt')).write_text(text)
        if name in ('mixture', 'signature') or smode:
            values[name] = '/work/localized/' + name + '.txt'
    env = WDL.values_from_json(values, task.available_inputs)
    stdlib = WDL.StdLib.Base('1.0')
    for decl in task.inputs:
        if not env.has_binding(decl.name):
            if decl.expr is not None:
                env = env.bind(decl.name, decl.expr.eval(env, stdlib))
            elif decl.type.optional:
                env = env.bind(decl.name, WDL.Value.Null())
    (base / 'command.sh').write_text(task.command.eval(env, stdlib).value)
    spec = importlib.util.spec_from_file_location('fraction_fixture', ROOT / 'tests/cibersortx_fractions/test_fractions.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    native = module.FractionsTest.native_fixture(0, None)
    native = native.replace("root = pathlib.Path(__file__).resolve().parent.parent", "root = pathlib.Path('/work')")
    (base / 'native_fixture.py').write_text(native)
    (base / 'native_fixture.py').chmod(0o755)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-only', type=Path)
    args = parser.parse_args()
    if args.prepare_only:
        args.prepare_only.mkdir(parents=True, exist_ok=True)
        prepare(args.prepare_only, True)
        print('Prepared S-mode command and native fixture; image not run.')
        return
    for smode in (False, True):
        with tempfile.TemporaryDirectory(prefix='fractions-image-') as temp:
            base = Path(temp)
            prepare(base, smode)
            # Check the real executable and Python runtime before replacing the
            # native calculation. This is an image compatibility check only.
            subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64',
                            '--entrypoint', '/bin/sh', IMAGE, '-c',
                            'test -x /src/CIBERSORTxFractions && python3 --version'], check=True)
            subprocess.run(['docker', 'run', '--rm', '--platform', 'linux/amd64',
                            '--mount', 'type=bind,src=%s,dst=/work' % base,
                            '--mount', 'type=bind,src=%s,dst=/src/CIBERSORTxFractions,readonly' % (base / 'native_fixture.py'),
                            '--workdir', '/work', '--entrypoint', '/bin/bash',
                            IMAGE, '-c',
                            # The image runs as root so it can stage /src paths.
                            # Return fixture files to the host runner before
                            # TemporaryDirectory removes them. Keep task errors.
                            'task_status=0; /bin/bash /work/command.sh || task_status=$?; '
                            'chown -R %d:%d /work || exit $?; exit "$task_status"' %
                            (os.getuid(), os.getgid())], check=True)
            report = json.loads((base / 'output_validation.json').read_text())
            assert report['sample_count'] == 3 and report['cell_types'] == ['B', 'CD4_T']
            assert report['smode'] == smode
            assert 'fixture-token' not in (base / 'fractions.log').read_text()
            print('PASS: image wrapper smoke test; S-mode=%s.' % smode)


if __name__ == '__main__':
    main()
