import csv
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
ENV = dict(os.environ, PYTHONPATH=str(ROOT / 'tools/trans_ld_regions/scripts'))

def table(path, fields, rows):
    with path.open('w', newline='') as h:
        w=csv.DictWriter(h, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(rows)

def read(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as h: return list(csv.DictReader(h, delimiter='\t'))

def run(module, args, cwd):
    return subprocess.run([sys.executable, '-m', 'trans_ld_regions.'+module, *map(str,args)], cwd=cwd, env=ENV, text=True, capture_output=True)

class PrepareTests(unittest.TestCase):
    def test_normalization_and_strict_input_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            a=p/'assoc.tsv'; a.write_text('variant_id\tphenotype_id\tpval\nchr1_100_A_G\tgeneA\t1e-9\nchr1_200_A_C\tgeneB\t0.1\n')
            (p/'files').write_text(str(a)+'\n')
            (p/'groups').write_text('IID\tancestry\nS1\tMID\nS2\tEUR\n')
            (p/'sizes').write_text('chr1\t10000000\n')
            (p/'chroms').write_text('1\n')
            args=['--association-list',p/'files','--ancestry',p/'groups','--chrom-sizes',p/'sizes','--chromosomes',p/'chroms','--pvalue-threshold','1e-5','--genome-build','GRCh38']
            r=run('prepare',args,p); self.assertEqual(r.returncode,0,r.stderr)
            rows=read(p/'associations.tsv.gz')
            self.assertEqual([(x['variant_id'],x['pos'],x['chrom']) for x in rows],[('chr1_100_A_G','100','1')])
            self.assertEqual((p/'ancestries.txt').read_text(),'EUR\nMID\n')
            (p/'groups').write_text('IID\tancestry\nS1\tMID\nS1\tEUR\n')
            self.assertNotEqual(run('prepare',args,p).returncode,0)
            (p/'files').write_text('gs://bucket/assoc.tsv\n')
            r=run('prepare',args,p); self.assertNotEqual(r.returncode,0); self.assertIn('localization error',r.stderr)

class MergeTests(unittest.TestCase):
    def test_union_padding_fallback_and_all_phenotype_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            assoc=[dict(variant_id=v,chrom='1',pos=str(pos),phenotype_id=ph,p_value='1e-9',modality='splicing',dataset_id='study',cell_type='CD4') for v,pos,ph in [('v1',500,'g1'),('v2',900,'g2'),('v3',5000,'g3')]]
            table(p/'associations.tsv',list(assoc[0]),assoc)
            intervals=[]
            for v,pos,group,lo,hi,status in [('v1',500,'EUR',300,600,'ok'),('v1',500,'MID',450,850,'ok'),('v2',900,'EUR',800,1000,'ok'),('v3',5000,'EUR',5000,5000,'missing_variant')]:
                intervals.append(dict(variant_id=v,chrom='1',pos=pos,ancestry=group,start=lo,end=hi,status=status,n_samples=100,search_bp=1000))
            table(p/'intervals.tsv',list(intervals[0]),intervals)
            (p/'files').write_text(str(p/'intervals.tsv')+'\n')
            (p/'sizes').write_text('1\t10000\n')
            r=run('combine',['--associations',p/'associations.tsv','--interval-list',p/'files','--chrom-sizes',p/'sizes','--padding-bp','100','--fallback-bp','1000','--min-region-bp','2000','--genome-build','GRCh38'],p)
            self.assertEqual(r.returncode,0,r.stderr)
            regions=read(p/'regions.tsv')
            self.assertEqual([(x['start'],x['end']) for x in regions],[('0','2000'),('3899','6100')])
            self.assertEqual(len(read(p/'region_associations.tsv.gz')),3)
            self.assertEqual(len(read(p/'trans_window_associations.tsv.gz')),3)
            self.assertIn('fallback',regions[1]['flags'])


class BoundaryTests(unittest.TestCase):
    def test_minimum_width_preserves_signal_and_shifts_at_chromosome_ends(self):
        sys.path.insert(0,str(ROOT/'tools/trans_ld_regions/scripts'))
        from trans_ld_regions.common import interval
        cases=[((3000000,3100000,10000000,2000000),(2050000,4050000)),
               ((0,100,10000000,2000000),(0,2000000)),
               ((9999900,10000000,10000000,2000000),(8000000,10000000)),
               ((1000000,5000000,10000000,2000000),(1000000,5000000)),
               ((10,100,1000,2000000),(0,1000))]
        for args,expected in cases:
            with self.subTest(args=args):self.assertEqual(interval(*args),expected)

    def test_empty_results_still_produce_valid_output_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            fields=['variant_id','chrom','pos','phenotype_id','p_value','modality','dataset_id','cell_type']
            table(p/'assoc.tsv',fields,[])
            (p/'files').write_text('')
            (p/'sizes').write_text('1\t10000000\n')
            r=run('combine',['--associations',p/'assoc.tsv','--interval-list',p/'files','--chrom-sizes',p/'sizes','--genome-build','GRCh38'],p)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual(read(p/'regions.tsv'),[])
            self.assertEqual(read(p/'region_associations.tsv.gz'),[])

class ParquetTests(unittest.TestCase):
    def test_parquet_and_tsv_inputs_preserve_distinct_source_context(self):
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError:
            self.skipTest('install pyarrow for Parquet input test')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            pq.write_table(pa.Table.from_pylist([dict(variant_id='rs1',chrom='chr1',pos=100,phenotype_id='g1',pval=1e-9,modality='expression')]),p/'input.parquet')
            (p/'input.tsv').write_text('variant_id\tchrom\tpos\tphenotype_id\tpval\nrs1\t1\t100\tg1\t1e-8\n')
            (p/'files').write_text(str(p/'input.parquet')+'\n'+str(p/'input.tsv')+'\n')
            (p/'groups').write_text('IID\tancestry\nS1\tEUR\n')
            (p/'sizes').write_text('chr1\t10000000\n');(p/'chroms').write_text('1\n')
            args=['--association-list',p/'files','--ancestry',p/'groups','--chrom-sizes',p/'sizes','--chromosomes',p/'chroms','--pvalue-threshold','1e-5','--genome-build','GRCh38']
            r=run('prepare',args,p);self.assertEqual(r.returncode,0,r.stderr)
            out=read(p/'associations.tsv.gz')
            self.assertEqual([r['dataset_id'] for r in out],['input_1','input_2'])
            self.assertEqual([r['pos'] for r in out],['100','100'])

if __name__ == '__main__':
    unittest.main()
