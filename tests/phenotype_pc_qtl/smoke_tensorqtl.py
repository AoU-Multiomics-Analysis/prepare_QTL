"""Run a small numerical smoke test inside the pinned upstream tensorQTL image.

No cloud jobs and no real participant data are used.
"""
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
import pgenlib

ROOT = Path(__file__).resolve().parents[2]
rng = np.random.RandomState(1234)
n = 80
samples = ['sample' + str(i) for i in range(n)]
genotypes = rng.binomial(2, .35, size=(4, n)).astype(np.int8)
with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp)
    # Localized File inputs can have unrelated paths and basenames.
    for name in ['pgen', 'pvar', 'psam']:
        (p / name).mkdir()
    pgen = p / 'pgen' / 'dosages.pgen'
    with pgenlib.PgenWriter(str(pgen).encode(), n, variant_ct=4) as writer:
        for row in genotypes:
            writer.append_biallelic(np.ascontiguousarray(row))
    pvar = p / 'pvar' / 'variants.pvar'
    pvar.write_text('#CHROM\tPOS\tID\tREF\tALT\nchr1\t100\tv1\tA\tC\nchrX\t100\tvX\tA\tC\nchrY\t100\tvY_near\tA\tC\nchrY\t8000000\tvY_far\tA\tC\n')
    psam = p / 'psam' / 'samples.psam'
    psam.write_text('#IID\n' + '\n'.join(samples) + '\n')
    cov = p / 'cov.tsv'
    values = pd.DataFrame([rng.normal(size=n), genotypes[0] + rng.normal(scale=.1, size=n), genotypes[1] + rng.normal(scale=.1, size=n)], index=['genPC1', 'exprPC1', 'exprPC2'], columns=samples)
    values.to_csv(cov, sep='\t', index_label='ID')
    ids = p / 'ids.txt'; ids.write_text('exprPC1\nexprPC2\n')
    for chromosome, expected in [('chrY', {'v1', 'vX', 'vY_far'}), ('PC_SCAN', {'v1', 'vX', 'vY_near', 'vY_far'})]:
        work = p / chromosome; work.mkdir()
        subprocess.run([sys.executable, str(ROOT / 'tools/phenotype_pc_qtl/prepare_pc_scan.py'), '--covariates', str(cov), '--phenotype-pc-ids', str(ids), '--chromosome', chromosome], cwd=work, check=True)
        subprocess.run([sys.executable, str(ROOT / 'tools/phenotype_pc_qtl/run_trans_pc.py'), '--plink-pgen', str(pgen), '--plink-pvar', str(pvar), '--plink-psam', str(psam), '--phenotype-bed', str(work / 'pc_scan.phenotypes.bed'), '--covariates', str(work / 'pc_scan.covariates.tsv'), '--maf-threshold', '0', '--pval-threshold', '1'], cwd=work, check=True)
        pairs = pd.read_parquet(work / 'pc_scan.trans_qtl_pairs.parquet')
        assert set(pairs.variant_id) == expected, pairs
        assert set(pairs.phenotype_id) == {'exprPC1', 'exprPC2'}
        assert len(pairs) == 2 * len(expected), pairs
        lead = pairs.loc[pairs.phenotype_id == 'exprPC1'].sort_values('pval').iloc[0]
        assert lead.variant_id == 'v1' and lead.pval < 1e-10, lead
        print('PASS: numerical tensorQTL scan and coordinate filtering:', chromosome, flush=True)
