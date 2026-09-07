"""CI only: render and execute each WDL task in the built image, without installing validators in it."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import WDL
from test_wdl import render
from test_ld import fixture

ROOT=Path(__file__).resolve().parents[2]

def main():
    image=sys.argv[1]
    doc=WDL.load(str(ROOT/'workflows/genotype/trans_ld_regions.wdl'))
    tasks={t.name:t for t in doc.tasks}
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp);p.chmod(0o777)
        data=p/"data ' $(touch INJECTED) `touch INJECTED2`";data.mkdir();fixture(data)
        def docker(command,cwd):
            subprocess.run(['docker','run','--rm','--user',str(os.getuid())+':'+str(os.getgid()),'-v',str(p)+':'+str(p),'-w',str(cwd),image,*command],check=True)
        docker(['plink2','--vcf',str(data/'input.vcf'),'--make-bed','--out',str(data/'geno')],p)
        paths={'gs://b/assoc':str(data/'assoc.tsv'),'gs://b/ancestry':str(data/'groups.tsv'),'gs://b/sizes':str(data/'sizes')}
        for ext in ['bed','bim','fam']:paths['gs://b/'+ext]=str(data/('geno.'+ext))
        def execute(name,values,d):
            d.mkdir()
            command=render(tasks[name],values,d,paths)
            (d/'command.sh').write_text(command)
            docker(['bash',str(d/'command.sh')],d)
        execute('PrepareInputs',dict(association_files=['gs://b/assoc'],ancestry_samples='gs://b/ancestry',chrom_sizes='gs://b/sizes',chromosome_names=['1'],genome_build='GRCh38',pvalue_threshold=1.0,pvalue_column='p_value',docker_image=image),p/'prepare')
        paths['gs://o/assoc']=str(p/'prepare'/'associations.tsv.gz')
        execute('AncestryLD',dict(genotype='gs://b/bed',variants='gs://b/bim',samples='gs://b/fam',genotype_format='bed',chromosome='1',ancestry_samples='gs://b/ancestry',ancestry_group='EUR',associations='gs://o/assoc',chrom_sizes='gs://b/sizes',docker_image=image,threads=1,memory_gb=2),p/'ld')
        paths['gs://o/intervals']=str(p/'ld'/'intervals.tsv')
        execute('CombineRegions',dict(associations='gs://o/assoc',interval_files=['gs://o/intervals'],chrom_sizes='gs://b/sizes',genome_build='GRCh38',docker_image=image),p/'merge')
        if not (p/'merge'/'regions.tsv').is_file():raise AssertionError('missing regions')
        if list(p.rglob('INJECTED')) or list(p.rglob('INJECTED2')):raise AssertionError('unsafe quoting')
        print('All three rendered WDL commands passed in the runtime image')

if __name__=='__main__':main()
