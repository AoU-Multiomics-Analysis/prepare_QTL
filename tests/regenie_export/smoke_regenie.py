"""Run real REGENIE steps on exports for all four molecular QTL types in CI."""
import csv
import gzip
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
RUNNER = '/opt/regenie/scripts/run_regenie.py'


def run(args, directory):
    subprocess.run([str(v) for v in args], cwd=directory, check=True)


def main():
    with tempfile.TemporaryDirectory(prefix='regenie_export_smoke_') as directory:
        root = Path(directory)
        rng = random.Random(7614)
        samples = [f's{i:04d}' for i in range(240)]
        dosage = []
        with (root / 'input.vcf').open('w') as out:
            out.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n')
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
            for variant in range(200):
                values = [int(rng.random() < .3) + int(rng.random() < .3) for _ in samples]
                dosage.append(values)
                out.write(f'{1+variant//100}\t{100+variant}\tv{variant}\tA\tG\t.\tPASS\t.\tGT\t' + '\t'.join(['0/0', '0/1', '1/1'][v] for v in values) + '\n')
        run(['plink2', '--vcf', root/'input.vcf', '--make-pgen', '--out', root/'geno', '--threads', '2'], root)
        # This matches the user's IID-only genotype inputs. The runner adds FID=0.
        (root/'geno.psam').write_text('#IID\n' + '\n'.join(samples) + '\n')
        (root/'common_variants.txt').write_text(''.join(f'v{i}\n' for i in range(200)))
        ages = [str(rng.uniform(20, 80)) for _ in samples]
        # Match the padded numeric text emitted by the existing R merge step.
        (root/'cov.tsv').write_text('ID\t' + '\t'.join(reversed(samples)) + '\nage\t' + '\t'.join(v.rjust(24) for v in reversed(ages)) + '\n')
        trait_values = [[str(.2*dosage[t*100][i]+rng.gauss(0,1)) for i in range(len(samples))] for t in range(2)]
        trait_values[0][1] = 'NA'
        for kind, ids in {'eqtl': ['ENSG001.1', 'ENSG002.2'], 'pqtl': ['protein:one', 'protein,two'],
                          'sqtl': ['chr1:2:3:clu_1', 'chr2:4:5:clu_2'], 'mqtl': ['chr1:100', 'chr2:200']}.items():
            work = root/kind
            work.mkdir()
            (work/'traits.bed').write_text('#chr\tstart\tend\tphenotype_id\t'+'\t'.join(samples)+'\n'+''.join(
                f'chr{t+1}\t0\t1\t{ids[t]}\t'+'\t'.join(trait_values[t])+'\n' for t in range(2)))
            run([sys.executable, ROOT/'tools/regenie_export/export_regenie.py', '--bed', work/'traits.bed',
                 '--covariates', root/'cov.tsv', '--covariate-format', 'tensor', '--output-prefix', kind], work)
            with (work/f'outputs/{kind}.batches.tsv').open() as stream:
                batch = next(csv.DictReader(stream, delimiter='\t'))
            common = ['--pgen', root/'geno.pgen', '--pvar', root/'geno.pvar', '--psam', root/'geno.psam',
                      '--phenotypes', work/'outputs'/batch['phenotype_file'], '--phenotype-columns', batch['phenotype_columns'],
                      '--covariates', work/f'outputs/{kind}.covariates.tsv', '--trait-type', 'quantitative',
                      '--threads', '2', '--block-size', '50']
            fit, test = work/'fit', work/'test'
            fit.mkdir(); test.mkdir()
            run([sys.executable, RUNNER, 'step1', *common, '--extract', root/'common_variants.txt',
                 '--lowmem', '--cv-folds', '3'], fit)
            predictions = sorted((fit/'results').glob('fit_*.loco'))
            (work/'prediction_files.list').write_text(''.join(str(p)+'\n' for p in predictions))
            run([sys.executable, RUNNER, 'step2', *common, '--prediction-list', fit/'results/fit_pred.list',
                 '--prediction-files', work/'prediction_files.list', '--min-mac', '5'], test)
            for trait in batch['phenotype_columns'].split(','):
                with gzip.open(test/f'results/assoc_{trait}.regenie.gz', 'rt') as stream:
                    result = [line.split() for line in stream]
                columns = result[0]
                p_index = columns.index('LOG10P')
                assert len(result) == 201, (kind, trait, len(result))
                assert any(row[p_index] != 'NA' and math.isfinite(float(row[p_index])) for row in result[1:])
            print(f'{kind}: exported phenotypes and covariates passed both REGENIE steps', flush=True)


if __name__ == '__main__':
    main()
