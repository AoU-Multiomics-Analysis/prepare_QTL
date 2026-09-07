"""Render actual WDL commands after simulated Cromwell File localization."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest
import WDL

ROOT=Path(__file__).resolve().parents[2]

class LocalStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename): return filename
    def _devirtualize_filename(self, filename): return filename


def render(task, values, directory, paths=None):
    env=WDL.values_from_json(values,task.available_inputs)
    if paths:
        env=WDL.Value.rewrite_env_paths(env,lambda f: paths.get(f.value,f.value))
    return task.command.eval(env,LocalStdLib('1.0',write_dir=str(directory/'writes'))).value

class WDLTests(unittest.TestCase):
    def test_localized_two_step_commands_and_quoting(self):
        doc=WDL.load(str(ROOT/'workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl'))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            # Deliberately include characters that would execute if quoting failed.
            unusual=p/"local ' $(touch INJECTED) `touch INJECTED2`"
            unusual.mkdir()
            cov=unusual/'cov.tsv'
            cov.write_text('ID\tS3\tS1\tS2\tS4\ngenPC1\t1\t2\t3\t4\nexprPC1\t4\t2\t1\t3\n')
            prep=p/'prepare'; prep.mkdir()
            command=render(doc.tasks[0],{'covariates':'gs://bucket/cov.tsv','phenotype_pc_ids':['exprPC1'],'placeholder_chromosome':'chrY','docker_image':'unused'},prep,{'gs://bucket/cov.tsv':str(cov)})
            r=subprocess.run(['bash','-c',command],cwd=prep,text=True,capture_output=True)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            task=doc.imports[0].doc.tasks[0]
            trans=p/'trans';trans.mkdir()
            paths={}
            values={'phenotype_bed':'gs://outputs/pc.bed','covariates':'gs://outputs/cov.tsv','maf_threshold':.05,'pval_threshold':1e-5,'num_threads':1,'memory_gb':2,'disk_gb':10,'preemptible_attempts':0,'docker_image':'unused'}
            paths['gs://outputs/pc.bed']=str(prep/'pc_scan.phenotypes.bed')
            paths['gs://outputs/cov.tsv']=str(prep/'pc_scan.covariates.tsv')
            contents={'pgen':'fixture','pvar':'#CHROM\tPOS\tID\tREF\tALT\n1\t100\tv1\tA\tC\n','psam':'#IID\nS1\nS2\nS3\nS4\n'}
            for ext,content in contents.items():
                d=unusual/ext;d.mkdir();f=d/('different_name.'+ext);f.write_text(content)
                uri='gs://genotypes/'+ext+'/file.'+ext
                values['plink_'+ext]=uri;paths[uri]=str(f)
            # Only tensorQTL computation is stubbed; actual embedded validation/staging code runs.
            module=p/'modules'/'tensorqtl';module.mkdir(parents=True)
            (module/'__init__.py').write_text('')
            (module/'__main__.py').write_text('from pathlib import Path\nimport sys\nbase=sys.argv[1]\nassert all(Path(base+"."+x).is_file() for x in ["pgen","pvar","psam"])\nassert "--mode" in sys.argv and "trans" in sys.argv\nassert Path(sys.argv[2]).is_file()\nassert Path(sys.argv[sys.argv.index("--covariates")+1]).is_file()\nPath("pc_scan.trans_qtl_pairs.parquet").write_text("stub")\n')
            command=render(task,values,trans,paths)
            env=dict(os.environ,PYTHONPATH=str(p/'modules'))
            r=subprocess.run(['bash','-c',command],cwd=trans,text=True,capture_output=True,env=env)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            self.assertTrue((trans/'pc_scan.trans_qtl_pairs.parquet').exists())
            self.assertFalse(list(p.rglob('INJECTED')))
            self.assertFalse(list(p.rglob('INJECTED2')))
            # A missed File localization must fail before tensorQTL runs.
            bad=p/'bad';bad.mkdir()
            r=subprocess.run(['bash','-c',render(task,values,bad)],cwd=bad,text=True,capture_output=True,env=env)
            self.assertNotEqual(r.returncode,0)
            self.assertIn('localization error',r.stdout+r.stderr)
            self.assertFalse((bad/'pc_scan.trans_qtl_pairs.parquet').exists())

    def test_generated_cloud_file_remains_file_typed(self):
        task=WDL.load(str(ROOT/'workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl')).tasks[0]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            cov=p/'cov.tsv'
            cov.write_text('ID\tA\tB\tC\tD\ngenPC1\t1\t2\t3\t4\nexprPC1\t4\t3\t1\t2\n')
            mapping={}
            class CloudStdLib(LocalStdLib):
                def _virtualize_filename(self, filename):
                    uri='gs://generated/' + Path(filename).name
                    mapping[uri]=filename
                    return uri
            env=WDL.values_from_json({'covariates':str(cov),'phenotype_pc_ids':['exprPC1'],'placeholder_chromosome':'chrY','docker_image':'unused'},task.available_inputs)
            stdlib=CloudStdLib('1.0',write_dir=str(p/'writes'))
            parts=[]
            for part in task.command.parts:
                if isinstance(part,str): parts.append(part)
                else:
                    value=part.expr.eval(env,stdlib)
                    # Model final command-time mapping of newly created File values.
                    if isinstance(value,WDL.Value.File):
                        value=WDL.Value.File(mapping.get(value.value,value.value))
                    parts.append(value.coerce(WDL.Type.String()).value)
            command=''.join(parts)
            self.assertNotIn('gs://generated/',command)
            run=subprocess.run(['bash','-c',command],cwd=p,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertTrue((p/'pc_scan.phenotypes.bed').is_file())

    def test_workflow_has_two_calls_and_no_file_writes(self):
        doc=WDL.load(str(ROOT/'workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl'))
        self.assertEqual(sum(isinstance(n,WDL.Tree.Call) for n in doc.workflow.body),2)
        self.assertEqual(doc.wdl_version,'1.0')
        def walk(node):
            if isinstance(node,WDL.Expr.Apply):
                self.assertNotIn(node.function_name,{'write_lines','write_tsv','write_map','write_json','write_objects','write_object'})
            for child in node.children: walk(child)
        walk(doc.workflow)
        for task in [doc.tasks[0],doc.imports[0].doc.tasks[0]]:
            file_names={'covariates'} if task.name=='PreparePCPhenotypes' else {'plink_pgen','plink_pvar','plink_psam','phenotype_bed','covariates'}
            for decl in task.inputs:
                if decl.name in file_names:self.assertIsInstance(decl.type,WDL.Type.File)

    def test_registration_and_external_image_policy(self):
        import yaml
        sys.path.insert(0,str(ROOT))
        from ci.plan_image_updates import plan_changes
        config=yaml.safe_load((ROOT/'ci/image-stages.yml').read_text())
        groups=[g for g in config['workflow_groups'] if g['paths']==['workflows/phenotype_pc_qtl/**']]
        self.assertEqual(len(groups),1)
        self.assertEqual(groups[0]['policy'],'external_images_manual')
        plan=plan_changes(config,['tools/phenotype_pc_qtl/prepare_pc_scan.py','workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl'])
        self.assertEqual(plan['builds'],[])
        self.assertEqual(plan['unmapped'],[])
        registration=yaml.safe_load((ROOT/'.dockstore.yml').read_text())
        entry=[x for x in registration['workflows'] if x['name']=='phenotype_pc_qtl']
        self.assertEqual(entry[0]['primaryDescriptorPath'],'/workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl')

    def test_embedded_scripts_match(self):
        subprocess.run([sys.executable,str(ROOT/'tools/phenotype_pc_qtl/sync_wdl.py'),'--check'],check=True)

if __name__=='__main__':unittest.main()
