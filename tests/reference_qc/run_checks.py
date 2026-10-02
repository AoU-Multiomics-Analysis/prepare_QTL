"""Check the published source snapshots without running scientific analysis."""
import ast
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'workflows/cell_type_specific_expression/reference_qc'


def main():
    record = json.loads((SOURCE / 'source_manifest.json').read_text())
    for entry in record['files']:
        path = ROOT / entry['published']
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError('Published source checksum changed: ' + str(path))
        if path.suffix == '.py':
            ast.parse(path.read_text(), filename=str(path))

    examples = ROOT / 'examples/reference_qc'
    doublets = json.loads((examples / 'scdblfinder.config.json').read_text())
    labels = json.loads((examples / 'singler.config.json').read_text())
    if (len(doublets['captures']) != 22 or
            doublets['captures'] != labels['captures'] or
            set(doublets['capture_seeds']) != set(doublets['captures'])):
        raise ValueError('Example capture lists and seeds must agree')
    if doublets['capture_seeds']['TSP7_Blood_NA_10X_1_1'] != 20261001:
        raise ValueError('The first TSP7 pilot seed changed')
    if doublets['capture_seeds']['TSP7_Blood_NA_10X_2_1'] != 20261002:
        raise ValueError('The second TSP7 pilot seed changed')

    for stage, test_dir in [('scdblfinder', '.'), ('singler', 'tests'),
                            ('singler_filter', 'tests')]:
        subprocess.run([sys.executable, '-m', 'unittest', 'discover',
                        '-s', test_dir, '-p', 'test_*.py', '-v'],
                       cwd=SOURCE / stage, check=True)
    print('Source integrity, example configurations, and 13 unit tests passed.', flush=True)


if __name__ == '__main__':
    main()
