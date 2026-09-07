"""Embed standard-library helpers so the WDL needs no script File inputs."""
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
for wdl, script in [('phenotype_pc_qtl.wdl', 'prepare_pc_scan.py'), ('tasks/tensorqtl_trans_compat.wdl', 'run_trans_pc.py')]:
    path = root / 'workflows/phenotype_pc_qtl' / wdl
    text = path.read_text()
    first, tail = text.split('# BEGIN EMBEDDED PYTHON\n')
    _, last = tail.split('# END EMBEDDED PYTHON\n')
    new = first + '# BEGIN EMBEDDED PYTHON\n' + (root / 'tools/phenotype_pc_qtl' / script).read_text() + '# END EMBEDDED PYTHON\n' + last
    if '--check' in sys.argv:
        if text != new:
            sys.exit('Embedded script differs: ' + wdl)
    else:
        path.write_text(new)
