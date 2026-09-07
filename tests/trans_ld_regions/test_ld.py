"""Real PLINK integration; enabled with PLINK2=/path/to/plink2."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from test_regions import run, read, table

PLINK = os.environ.get('PLINK2') or shutil.which('plink2')


def fixture(p):
    ids = [g+str(i) for g in ['EUR','MID'] for i in range(100)]
    (p/'groups.tsv').write_text('IID\tancestry\n'+''.join(s+'\t'+('EUR' if s.startswith('EUR') else 'MID')+'\n' for s in ids))
    with (p/'input.vcf').open('w') as h:
        h.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=20000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        h.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n')
        for pos,v in [(3900,'rare'),(4100,'left'),(5000,'seed'),(5900,'bridge'),(6900,'right'),(8900,'far')]:
            gt=[]
            for group in ['EUR','MID']:
                for i in range(100):
                    a=['0/0','0/1','1/1','0/1'][i%4]
                    b=['0/0','0/1','1/1','0/1'][(i//4)%4]
                    if v=='rare': val='0/1' if i<2 else '0/0'
                    elif v=='seed': val=a
                    elif v=='left': val=a if group=='EUR' else b
                    else: val=a if group=='MID' else b
                    gt.append(val)
            h.write(f'1\t{pos}\t{v}\tA\tG\t.\tPASS\t.\tGT\t'+'\t'.join(gt)+'\n')
    assoc=[dict(variant_id=v,chrom='1',pos=pos,phenotype_id=v,p_value='1e-9',modality='splicing',dataset_id='test',cell_type='test') for v,pos in [('seed',5000),('rare',3900),('missing',15000)]]
    table(p/'assoc.tsv',list(assoc[0]),assoc)
    (p/'sizes').write_text('1\t20000\n')
    return ids

@unittest.skipUnless(PLINK,'set PLINK2 for real PLINK smoke tests')
class LDTests(unittest.TestCase):
    def test_bed_pgen_ancestry_ld_strict_maf_and_search_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); fixture(p)
            for mode,exts in [('bed',['bed','bim','fam']),('pgen',['pgen','pvar','psam'])]:
                command=[PLINK,'--vcf',str(p/'input.vcf'),'--make-'+mode,'--out',str(p/mode)]
                r=subprocess.run(command,text=True,capture_output=True)
                self.assertEqual(r.returncode,0,r.stdout+r.stderr)
                for group in ['EUR','MID']:
                    d=p/(mode+group); d.mkdir()
                    args=['--genotype',p/(mode+'.'+exts[0]),'--variants',p/(mode+'.'+exts[1]),'--samples',p/(mode+'.'+exts[2]),'--format',mode,'--chromosome','1','--ancestry',p/'groups.tsv','--group',group,'--associations',p/'assoc.tsv','--chrom-sizes',p/'sizes','--r2-threshold','0.8','--initial-search-bp','1000','--max-search-bp','4000','--edge-bp','100','--plink',PLINK,'--threads','1','--memory-mb','1024']
                    result=run('ld',args,d)
                    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                    out={x['variant_id']:x for x in read(d/'intervals.tsv')}
                    self.assertEqual(out['rare']['status'],'maf_filtered')
                    self.assertEqual(float(out['rare']['seed_mac']),2)
                    self.assertEqual(out['missing']['status'],'missing_variant')
                    self.assertEqual(int(out['seed']['n_samples']),100)
                    if group=='EUR':
                        self.assertEqual((out['seed']['start'],out['seed']['end']),('4100','5000'))
                        self.assertEqual(out['seed']['status'],'ok')
                    else:
                        self.assertEqual(out['seed']['status'],'search_limit')
                        self.assertEqual((out['seed']['start'],out['seed']['end']),('1000','9000'))
                        self.assertEqual(out['seed']['right_variant'],'far')


@unittest.skipUnless(PLINK,'set PLINK2 for real PLINK smoke tests')
class EdgeLDTests(unittest.TestCase):
    def test_numeric_x_code_and_constant_dosage(self):
        for constant in [False,True]:
            with self.subTest(constant=constant), tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp);fixture(p)
                if constant:
                    text=(p/'input.vcf').read_text().splitlines()
                    text=['\t'.join(line.split('\t')[:9]+['0/1']*200) if not line.startswith('#') and line.split('\t')[2]=='seed' else line for line in text]
                    (p/'input.vcf').write_text('\n'.join(text)+'\n')
                r=subprocess.run([PLINK,'--vcf',str(p/'input.vcf'),'--make-bed','--out',str(p/'geno')],capture_output=True,text=True)
                self.assertEqual(r.returncode,0,r.stderr)
                # Standard PLINK 1 encoding: chromosome X is numeric 23 in BIM.
                (p/'geno.bim').write_text(''.join('23\t'+'\t'.join(line.split()[1:])+'\n' for line in (p/'geno.bim').read_text().splitlines()))
                fam=[]
                for line in (p/'geno.fam').read_text().splitlines():
                    f=line.split();f[4]='2';fam.append('\t'.join(f))
                (p/'geno.fam').write_text('\n'.join(fam)+'\n')
                a=read(p/'assoc.tsv')
                for row in a:row['chrom']='X'
                table(p/'assoc.tsv',list(a[0]),a)
                (p/'sizes').write_text('X\t20000\n')
                args=['--genotype',p/'geno.bed','--variants',p/'geno.bim','--samples',p/'geno.fam','--format','bed','--chromosome','X','--ancestry',p/'groups.tsv','--group','EUR','--associations',p/'assoc.tsv','--chrom-sizes',p/'sizes','--plink',PLINK,'--threads','1','--memory-mb','1024']
                r=run('ld',args,p);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
                out={x['variant_id']:x for x in read(p/'intervals.tsv')}
                self.assertEqual(out['seed']['status'],'no_reported_ld' if constant else 'ok')

@unittest.skipUnless(PLINK,'set PLINK2 for real PLINK smoke tests')
class SparseLDTests(unittest.TestCase):
    def test_low_mac_is_retained_and_small_group_is_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);fixture(p)
            text=(p/'input.vcf').read_text().splitlines()
            # 3/200 allele copies = 1.5% in each group. A MAC20 filter would remove it.
            for i,line in enumerate(text):
                if not line.startswith('#') and line.split('\t')[2]=='rare':
                    text[i]='\t'.join(line.split('\t')[:9]+[('0/1' if j%100<3 else '0/0') for j in range(200)])
            (p/'input.vcf').write_text('\n'.join(text)+'\n')
            r=subprocess.run([PLINK,'--vcf',str(p/'input.vcf'),'--make-bed','--out',str(p/'geno')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            args=['--genotype',p/'geno.bed','--variants',p/'geno.bim','--samples',p/'geno.fam','--format','bed','--chromosome','1','--ancestry',p/'groups.tsv','--group','EUR','--associations',p/'assoc.tsv','--chrom-sizes',p/'sizes','--plink',PLINK,'--threads','1','--memory-mb','1024']
            d=p/'full';d.mkdir();r=run('ld',args,d);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            freq={x['variant_id']:x for x in read(d/'allele_counts.tsv.gz')}
            self.assertEqual(freq['rare']['eligible'],'True');self.assertEqual(float(freq['rare']['mac']),3)
            (p/'groups.tsv').write_text('IID\tancestry\n'+''.join('EUR'+str(i)+'\tEUR\n' for i in range(20)))
            d=p/'small';d.mkdir();r=run('ld',args,d);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            out={x['variant_id']:x for x in read(d/'intervals.tsv')}
            self.assertEqual(out['seed']['status'],'insufficient_samples')
            self.assertEqual(out['seed']['n_samples'],'20')

@unittest.skipUnless(PLINK,'set PLINK2 for real PLINK smoke tests')
class CompressedPvarTests(unittest.TestCase):
    def test_compressed_pvar_and_mismatched_coordinates(self):
        try:
            import zstandard
        except ImportError:
            self.skipTest('install zstandard for compressed PVAR test')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);fixture(p)
            r=subprocess.run([PLINK,'--vcf',str(p/'input.vcf'),'--make-pgen','--out',str(p/'geno')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            (p/'geno.pvar.zst').write_bytes(zstandard.ZstdCompressor().compress((p/'geno.pvar').read_bytes()))
            args=['--genotype',p/'geno.pgen','--variants',p/'geno.pvar.zst','--samples',p/'geno.psam','--format','pgen','--chromosome','1','--ancestry',p/'groups.tsv','--group','EUR','--associations',p/'assoc.tsv','--chrom-sizes',p/'sizes','--plink',PLINK,'--threads','1','--memory-mb','1024']
            d=p/'good';d.mkdir();r=run('ld',args,d);self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            self.assertEqual(next(x for x in read(d/'intervals.tsv') if x['variant_id']=='seed')['status'],'ok')
            a=read(p/'assoc.tsv');a[0]['pos']='5001';table(p/'assoc.tsv',list(a[0]),a)
            d=p/'bad';d.mkdir();r=run('ld',args,d);self.assertNotEqual(r.returncode,0)
            self.assertIn('position mismatch',r.stderr)

if __name__ == '__main__':
    unittest.main()
