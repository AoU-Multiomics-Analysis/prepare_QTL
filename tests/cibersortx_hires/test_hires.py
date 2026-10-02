"""Test the WDL command at its file and native-process boundaries.

Cloud downloads are simulated by replacing File values with readable local
paths, as a workflow engine must do. This is not a Terra integration test.
Only the authenticated native calculation is replaced with a small executable.
"""
import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

ROOT = Path(__file__).resolve().parents[2]
WDL_PATH = ROOT / "workflows/cell_type_specific_expression/cibersortx_hires.wdl"


class HiResTest(unittest.TestCase):
    def document(self):
        path = WDL_PATH
        self.assertTrue(path.is_file(), "The Terra WDL has not been created")
        return WDL.load(str(path))

    def run_task(self, change=None, native_status=0, bad_output=False, inspect=None):
        doc = self.document()
        task = doc.tasks[0]
        with tempfile.TemporaryDirectory(prefix="hires localized ") as tmp:
            base = Path(tmp)
            localized = base / "localized input's"
            localized.mkdir()
            contents = {
                "mixture": "GeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nsig\t1\t2\t3\t4\t5\t6\t7\t8\nZNF804A\t2\t3\t4\t5\t6\t7\t8\t9\n",
                "signature": "GeneSymbol\tB\tCD4_T\nsig\t1\t2\n",
                "fractions": "Mixture\tB\tCD4_T\ns1\t0.2\t0.8\ns2\t0.5\t0.5\ns3\t0.8\t0.2\ns4\t0.3\t0.7\ns5\t0.4\t0.6\ns6\t0.6\t0.4\ns7\t0.7\t0.3\ns8\t0.9\t0.1\n",
                "gene_subset": "ZNF804A\n",
            }
            inputs = {}
            for name, content in contents.items():
                file = localized / (name + ".txt")
                file.write_text(content)
                inputs[name] = str(file)
            inputs.update(username="user'$(touch INJECTED)@example.org", token="token'`touch INJECTED`", threads=8,
                          nsampling=1, nsampling2=1, memory_gb=16, disk_gb=20)
            if change:
                change(inputs, localized)
            # Start with cloud File inputs. Simulate the engine's download by
            # rewriting File values, while leaving String values untouched.
            cloud_inputs = dict(inputs)
            downloads = {}
            for name in contents:
                if '://' not in inputs[name]:
                    cloud_uri = 'gs://test-bucket/' + name + '.txt'
                    downloads[cloud_uri] = inputs[name]
                    cloud_inputs[name] = cloud_uri
            for name in contents:
                self.assertIsInstance(doc.workflow.available_inputs.resolve(name).type, WDL.Type.File)
            env = WDL.values_from_json(cloud_inputs, task.available_inputs)
            env = WDL.Value.rewrite_env_paths(env, lambda file: downloads.get(file.value, file.value))
            stdlib = WDL.StdLib.Base("1.0")
            for decl in task.inputs:
                if not env.has_binding(decl.name) and decl.expr is not None:
                    env = env.bind(decl.name, decl.expr.eval(env, stdlib))
            runtime_cpu = task.runtime["cpu"].eval(env, stdlib).value
            command = task.command.eval(env, stdlib).value
            # Preserve the command. Map only the image's fixed /src directory
            # into the temporary test filesystem.
            src = base / "src"
            src.mkdir()
            native = src / "CIBERSORTxHiRes"
            native.write_text("#!/usr/bin/env python3\n" +
                "import json, pathlib, sys\n" +
                "args = dict(zip(sys.argv[1::2], sys.argv[2::2]))\n" +
                "out = pathlib.Path(args['--outdir']).resolve()\n" +
                "(out.parent / 'native_args.json').write_text(json.dumps(args))\n" +
                "print('Native calculation started', flush=True)\n" +
                ("sys.exit(%d)\n" % native_status if native_status else
                 "for cell in ['B', 'CD4_T']:\n" +
                 "    (out / ('CIBERSORTxHiRes_' + args['--label'] + '_' + cell + '_Window8.txt')).write_text(" +
                 repr("GeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\ts7\t" + ("WRONG" if bad_output else "s8") + "\nZNF804A\t2\t3\t4\t5\t6\t7\t8\t9\n") + ")\n"))
            native.chmod(0o755)
            image = os.environ.get('CIBERSORTX_TEST_IMAGE')
            # Local tests map /src into a temporary directory. Container tests
            # use the actual image paths. Input File values stay intact.
            if not image:
                command = command.replace("/src", str(src))
            script = base / "command.sh"
            script.write_text(command)
            invocation = ["bash", str(script)]
            if image:
                invocation = ['docker', 'run', '--rm', '--platform', 'linux/amd64',
                              '--entrypoint', '/bin/bash', '--workdir', str(base),
                              '-v', str(base) + ':' + str(base),
                              '-v', str(native) + ':/src/CIBERSORTxHiRes:ro', image,
                              '-c', 'set +e; bash "$1"; result=$?; chown -R "$2:$3" "$4"; exit "$result"',
                              'hires-test', str(script), str(os.getuid()), str(os.getgid()), str(base)]
            proc = subprocess.run(invocation, cwd=base, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
            args_path = base / "native_args.json"
            args = json.loads(args_path.read_text()) if args_path.exists() else None
            report = base / "output_validation.json"
            if inspect:
                inspect(base)
            return proc, args, runtime_cpu, bool(list(base.rglob("INJECTED"))), json.loads(report.read_text()) if report.exists() else None

    def test_localized_files_and_literal_flags_reach_native_process(self):
        # Break caught: quoting loses apostrophes, executes shell text, or
        # typed File values are replaced with unresolved cloud metadata.
        proc, args, cpu, injected, report = self.run_task()
        self.assertEqual(proc.returncode, 0, proc.stdout)
        self.assertEqual(args["--username"], "user'$(touch INJECTED)@example.org")
        self.assertEqual(args["--token"], "token'`touch INJECTED`")
        self.assertEqual((cpu, args["--threads"]), (8, "8"))
        self.assertFalse(injected)
        self.assertEqual((report["sample_count"], report["cell_types"]), (8, ["B", "CD4_T"]))

    def test_unresolved_cloud_file_fails_before_native_process(self):
        for name in ("mixture", "fractions"):
            with self.subTest(file=name):
                proc, args, _, _, _ = self.run_task(
                    lambda inputs, _: inputs.update({name: "gs://bucket/" + name + ".txt"}))
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("Localization error: " + name, proc.stdout)
                self.assertIsNone(args)

    def test_multiline_flag_is_passed_as_one_literal_argument(self):
        # Break caught: a rendered String terminates a shell wrapper here-doc.
        token = "first'line\nHIRES_TASK\nlast`touch INJECTED`"
        proc, args, _, injected, _ = self.run_task(lambda inputs, _: inputs.update(token=token))
        self.assertEqual(proc.returncode, 0, proc.stdout)
        self.assertEqual(args['--token'], token)
        self.assertFalse(injected)

    def test_fit_metric_columns_fail_before_native_process(self):
        def metrics(inputs, _):
            Path(inputs["fractions"]).write_text("Mixture\tB\tCD4_T\tP-value\tCorrelation\tRMSE\ns1\t0.2\t0.8\t0\t0.9\t0.1\n")
        proc, args, _, _, _ = self.run_task(metrics)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("fraction-only", proc.stdout)
        self.assertIsNone(args)

    def test_fraction_samples_are_reordered_to_match_mixture(self):
        # Break caught: reversed rows are rejected or values stay in donor order.
        expected = (
            "Mixture\tB\tCD4_T\n"
            "s1\t0.2\t0.8\ns2\t0.5\t0.5\ns3\t0.8\t0.2\ns4\t0.3\t0.7\n"
            "s5\t0.4\t0.6\ns6\t0.6\t0.4\ns7\t0.7\t0.3\ns8\t0.9\t0.1\n"
        )
        original = {}
        observed = {}
        def reorder(inputs, _):
            file = Path(inputs['fractions'])
            rows = file.read_text().splitlines()
            file.write_text('\n'.join([rows[0]] + rows[:0:-1]) + '\n')
            original['bytes'] = file.read_bytes()
        def inspect(base):
            observed['original'] = (base / "localized input's/fractions.txt").read_bytes()
            staged = base / 'results/fractions.txt'
            observed['staged'] = staged.read_text() if staged.exists() else None
            report = base / 'input_validation.json'
            observed['report'] = json.loads(report.read_text()) if report.exists() else {}
        proc, args, _, _, _ = self.run_task(reorder, inspect=inspect)
        self.assertEqual(proc.returncode, 0, proc.stdout)
        self.assertIsNotNone(args)
        self.assertEqual(args['--cibresults'], 'fractions.txt')
        self.assertEqual(observed['staged'], expected)
        self.assertEqual(observed['original'], original['bytes'])
        self.assertEqual(observed['report']['fraction_input_sample_count'], 8)
        self.assertEqual(observed['report']['fraction_retained_sample_count'], 8)
        self.assertEqual(observed['report']['fraction_dropped_sample_count'], 0)

    def test_fraction_superset_is_selected_without_changing_numeric_text(self):
        # Break caught: extra donors enter HiRes, or float serialization changes values.
        source = (
            "Sample ID\tB\tCD4_T\n"
            "extra_a\t0.25\t0.75\ns8\t9e-1\t1e-1\ns7\t.700\t.300\n"
            "s6\t0.6000\t0.4000\nextra_b\t0.75\t0.25\ns5\t4e-1\t6e-1\n"
            "s4\t.30\t.70\ns3\t0.80000\t0.20000\ns2\t5.0e-1\t5.0e-1\n"
            "extra_c\t0.1\t0.9\ns1\t0.200000\t0.800000\n"
        )
        expected = (
            "Sample ID\tB\tCD4_T\n"
            "s1\t0.200000\t0.800000\ns2\t5.0e-1\t5.0e-1\ns3\t0.80000\t0.20000\n"
            "s4\t.30\t.70\ns5\t4e-1\t6e-1\ns6\t0.6000\t0.4000\n"
            "s7\t.700\t.300\ns8\t9e-1\t1e-1\n"
        )
        observed = {}
        def superset(inputs, _):
            Path(inputs['fractions']).write_text(source)
        def inspect(base):
            observed['original'] = (base / "localized input's/fractions.txt").read_bytes()
            staged = base / 'results/fractions.txt'
            observed['staged'] = staged.read_text() if staged.exists() else None
            report = base / 'input_validation.json'
            observed['report'] = json.loads(report.read_text()) if report.exists() else {}
        proc, args, _, _, _ = self.run_task(superset, inspect=inspect)
        self.assertEqual(proc.returncode, 0, proc.stdout)
        self.assertIsNotNone(args)
        self.assertEqual(args['--cibresults'], 'fractions.txt')
        self.assertEqual(observed['staged'], expected)
        self.assertEqual(observed['original'], source.encode('utf-8'))
        self.assertEqual(observed['report']['fraction_input_sample_count'], 11)
        self.assertEqual(observed['report']['fraction_retained_sample_count'], 8)
        self.assertEqual(observed['report']['fraction_dropped_sample_count'], 3)
        self.assertEqual(observed['report']['samples'], ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8'])
        self.assertRegex(proc.stdout, r'(?i)(?:selected|retained)\D*8\b')
        self.assertRegex(proc.stdout, r'(?i)dropped\D*3\b')

    def test_missing_fraction_sample_fails_before_native_process(self):
        # Break caught: a missing donor is silently removed or replaced by an extra row.
        def missing(inputs, _):
            file = Path(inputs['fractions'])
            rows = file.read_text().splitlines()
            file.write_text('\n'.join(rows[:-1] + ['extra\t0.5\t0.5']) + '\n')
        proc, args, _, _, _ = self.run_task(missing)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn('s8', proc.stdout)
        self.assertIsNone(args)

    def test_duplicate_fraction_samples_fail_before_native_process(self):
        # Break caught: indexing rows silently overwrites retained or excluded duplicates.
        for duplicate in ('s1\t0.2\t0.8\n', 'extra\t0.2\t0.8\nextra\t0.8\t0.2\n'):
            with self.subTest(rows=duplicate):
                def add_duplicate(inputs, _):
                    file = Path(inputs['fractions'])
                    file.write_text(file.read_text() + duplicate)
                proc, args, _, _, _ = self.run_task(add_duplicate)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn('duplicate labels', proc.stdout)
                self.assertIsNone(args)

    def test_cell_type_with_zero_retained_fraction_fails_before_native_process(self):
        # Break caught: an excluded donor supplies the only positive cell-type fraction.
        def zero_retained(inputs, _):
            Path(inputs['fractions']).write_text(
                'Mixture\tB\tCD4_T\nextra\t0\t1\n' +
                ''.join('s%d\t1\t0\n' % sample for sample in range(1, 9)))
        proc, args, _, _, _ = self.run_task(zero_retained)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn('positive total fraction', proc.stdout)
        self.assertIsNone(args)

    def test_native_failure_is_not_hidden_by_logging(self):
        proc, _, _, _, _ = self.run_task(native_status=17)
        self.assertEqual(proc.returncode, 17, proc.stdout)
        self.assertIn("Native calculation started", proc.stdout)

    def test_wrong_output_samples_fail_the_task(self):
        proc, _, _, _, _ = self.run_task(bad_output=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("Output error: sample columns", proc.stdout)

    def test_workflow_file_writer_check(self):
        path = ROOT / "scripts/check_wdl_file_scope.py"
        self.assertTrue(path.is_file(), "The workflow-scope regression check has not been created")
        spec = importlib.util.spec_from_file_location("check_workflow", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        forbidden = WDL.parse_document('version 1.0\nworkflow bad { Array[String] x = ["x"] File f = write_lines(x) }')
        self.assertTrue(module.workflow_file_writes(forbidden.workflow))
        nested = WDL.parse_document('version 1.0\ntask consumer { input { String text } command { echo ~{text} } }\nworkflow bad { call consumer { input: text=read_string(write_lines(["x"])) } }')
        self.assertTrue(module.workflow_file_writes(nested.workflow))
        task_only = WDL.parse_document('version 1.0\ntask ok { Array[String] x = ["x"] File f = write_lines(x) command { cat ~{f} } }\nworkflow good { call ok }')
        self.assertEqual(module.workflow_file_writes(task_only.workflow), [])
        self.assertEqual(module.workflow_file_writes(self.document().workflow), [])


if __name__ == "__main__":
    unittest.main()
