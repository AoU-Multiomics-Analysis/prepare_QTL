"""Terra-oriented WDL checks and real rendered commands after File localization."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import WDL
from test_regions import ROOT, read
from test_ld import fixture, PLINK

class LocalStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename): return filename
    def _devirtualize_filename(self, filename): return filename

def render(task, values, directory, paths):
    env=WDL.values_from_json(values,task.available_inputs)
    env=WDL.Value.rewrite_env_paths(env,lambda f:paths.get(f.value,f.value))
    lib=LocalStdLib('1.0',write_dir=str(directory/'lists'))
    for decl in task.inputs:
        if not env.has_binding(decl.name) and decl.expr:
            env=env.bind(decl.name,decl.expr.eval(env,lib))
    return task.command.eval(env,lib).value

class WDLTests(unittest.TestCase):
    def test_no_workflow_scope_file_writes_and_file_inputs(self):
        self.assertTrue((ROOT/'workflows/genotype/trans_ld_regions.wdl').exists(),'workflow has not been implemented')
        doc=WDL.load(str(ROOT/'workflows/genotype/trans_ld_regions.wdl'))
        self.assertEqual(doc.wdl_version,'1.0')
        def walk(n):
            if isinstance(n,WDL.Expr.Apply):
                self.assertNotIn(n.function_name,{'write_lines','write_tsv','write_map','write_json','write_object','write_objects'})
            for child in n.children: walk(child)
        walk(doc.workflow)
        for task in doc.tasks:
            for d in task.inputs:
                if d.name in {'ancestry_samples','associations','chrom_sizes','genotype','variants','samples'}:
                    self.assertIsInstance(d.type,WDL.Type.File)
                if d.name in {'association_files','interval_files'}:
                    self.assertIsInstance(d.type.item_type,WDL.Type.File)

    @unittest.skipUnless(PLINK,'set PLINK2 for rendered-command integration')
    def test_localization_all_tasks_and_safe_quoting(self):
        self.assertTrue((ROOT/'workflows/genotype/trans_ld_regions.wdl').exists(),'workflow has not been implemented')
        doc=WDL.load(str(ROOT/'workflows/genotype/trans_ld_regions.wdl'))
        tasks={t.name:t for t in doc.tasks}
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); data=p/"data ' $(touch INJECTED) `touch INJECTED2`"; data.mkdir(); fixture(data)
            r=subprocess.run([PLINK,'--vcf',str(data/'input.vcf'),'--make-bed','--out',str(data/'geno')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            paths={'gs://b/assoc.tsv':str(data/'assoc.tsv'),'gs://b/groups.tsv':str(data/'groups.tsv'),'gs://b/sizes':str(data/'sizes')}
            for ext in ['bed','bim','fam']: paths['gs://g/file.'+ext]=str(data/('geno.'+ext))
            env=dict(os.environ,PYTHONPATH=str(ROOT/'tools/trans_ld_regions/scripts'),PATH=str(Path(PLINK).parent)+os.pathsep+os.environ['PATH'])
            prep=p/'prep'; prep.mkdir()
            pv=dict(association_files=['gs://b/assoc.tsv'],ancestry_samples='gs://b/groups.tsv',chrom_sizes='gs://b/sizes',chromosome_names=['1'],pvalue_threshold=1e-5,genome_build='GRCh38',pvalue_column='p_value',docker_image='unused')
            def execute(name,values,d,localized=True):
                result=subprocess.run(['bash','-c',render(tasks[name],values,d,paths if localized else {})],cwd=d,env=env,text=True,capture_output=True)
                return result
            r=execute('PrepareInputs',pv,prep); self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            paths['gs://o/associations']=str(prep/'associations.tsv.gz')
            ld=p/'ld'; ld.mkdir()
            lv=dict(genotype='gs://g/file.bed',variants='gs://g/file.bim',samples='gs://g/file.fam',genotype_format='bed',chromosome='1',ancestry_samples='gs://b/groups.tsv',ancestry_group='EUR',associations='gs://o/associations',chrom_sizes='gs://b/sizes',docker_image='unused',threads=1,memory_gb=2)
            r=execute('AncestryLD',lv,ld); self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            paths['gs://o/intervals']=str(ld/'intervals.tsv')
            merge=p/'merge';merge.mkdir()
            mv=dict(associations='gs://o/associations',interval_files=['gs://o/intervals'],chrom_sizes='gs://b/sizes',genome_build='GRCh38',docker_image='unused')
            r=execute('CombineRegions',mv,merge); self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            self.assertEqual(len(read(merge/'region_associations.tsv.gz')),3)
            for name,values in [('PrepareInputs',pv),('AncestryLD',lv),('CombineRegions',mv)]:
                d=p/('bad_'+name);d.mkdir()
                r=execute(name,values,d,False)
                self.assertNotEqual(r.returncode,0); self.assertIn('localization error',r.stdout+r.stderr)
            self.assertFalse(list(p.rglob('INJECTED')))
            self.assertFalse(list(p.rglob('INJECTED2')))


class GeneratedFileTests(unittest.TestCase):
    def test_new_cloud_file_results_remain_files_and_can_be_localized(self):
        class CloudStdLib(LocalStdLib):
            def __init__(self,directory):
                super().__init__('1.0',write_dir=str(directory));self.generated={}
            def _virtualize_filename(self,filename):
                uri='gs://generated/'+Path(filename).name;self.generated[uri]=filename;return uri
        doc=WDL.load(str(ROOT/'workflows/genotype/trans_ld_regions.wdl'))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            values={'association_files':['gs://input/assoc'], 'chromosome_names':['1'], 'interval_files':['gs://input/interval']}
            for task in doc.tasks:
                for part in task.command.parts:
                    if isinstance(part,WDL.Expr.Placeholder) and 'write_lines' in str(part.expr):
                        names={d.name for d in task.inputs}
                        env=WDL.values_from_json({k:v for k,v in values.items() if k in names},task.available_inputs)
                        env=WDL.Value.rewrite_env_paths(env,lambda f:'/localized/'+Path(f.value).name)
                        lib=CloudStdLib(p)
                        value=part.expr.eval(env,lib)
                        self.assertIsInstance(value,WDL.Value.File,'generated list lost its File type')
                        local=WDL.Value.rewrite_paths(value,lambda f:lib.generated[f.value])
                        expected=['1'] if 'chromosome_names' in str(part.expr) else ['/localized/assoc' if 'association_files' in str(part.expr) else '/localized/interval']
                        self.assertEqual(Path(local.value).read_text().splitlines(),expected)

if __name__ == '__main__':
    unittest.main()
