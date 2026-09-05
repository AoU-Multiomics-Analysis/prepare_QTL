"""Static routing checks use candidate WDL and independent trusted policy."""
from pathlib import Path
import sys
import tempfile
import unittest

import WDL
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ci'))
from wdl_stage_routing import task_stages, validate_routing


class WdlStageRoutingTest(unittest.TestCase):
    TRANSFER_PROGRAMS = (
        'gsutil cp "$1" "$2"',
        '''
        sample_id="$1"
        source_path="$2"
        local_path="$3"
        if [[ "$source_path" == gs://* ]]; then
            gsutil -q cp "$source_path" "$local_path"
        else
            if [ ! -f "$source_path" ]; then
                echo "Input BED file for ${sample_id} is not accessible inside the task: ${source_path}" >&2
                exit 1
            fi
            cp "$source_path" "$local_path"
        fi
        ''',
    )

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'source'
        self.policy = Path(self.temp.name) / 'policy'
        self.config = {
            'images': {'cell': {'repository': 'ghcr.io/example/cell'}},
            'stages': {
                'cell_export': {'image': 'cell', 'script_roots': ['scripts/cell_type_specific_expression/export']},
                'cell_fit': {'image': 'cell', 'script_roots': ['scripts/cell_type_specific_expression/fit']},
            },
            'shared': [{'sources': ['scripts/common/**'], 'stages': ['cell_export', 'cell_fit']}],
        }
        self.pins = {'stages': {s: [] for s in self.config['stages']}}
        self.script('export/new.R')
        self.script('fit/fit.R')
        self.write(self.root, 'scripts/common/helper.R', '# test only')
        self.command = 'Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/export/new.R'

    def write(self, root, name, text):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def script(self, path):
        self.write(self.root, 'scripts/cell_type_specific_expression/' + path, '# test only')

    def workflow(self, path='workflows/main.wdl', binding='export_docker_image',
                 child=False, omit=False, runtime='docker_image', default=False,
                 wrapper='', command=None):
        inputs = '\n'.join(f'String {s}_docker_image = "same-image"' for s in ('export', 'fit'))
        for short in ('export', 'fit'):
            self.pins['stages']['cell_' + short].append({'path': path, 'input': short + '_docker_image'})
        if child:
            call = 'call child.Main as RunChild { input: export_docker_image = ' + binding + ', fit_docker_image = fit_docker_image }'
            header = 'import "child.wdl" as child\n'
            task = ''
        else:
            call = 'call NeverRegistered' + ('' if omit else ' { input: docker_image = ' + binding + ' }')
            header = ''
            task = '''task NeverRegistered {
 input { String docker_image%s String script = "x" }
 command <<< %s >>>
 runtime { docker: %s }
}''' % (' = "same-image"' if default else '', command or self.command, runtime)
        if wrapper:
            call = wrapper + ' { ' + call + ' }'
        self.write(self.root, path, 'version 1.0\n' + header + task + '\nworkflow Main { input { ' + inputs + ' } ' + call + ' }')

    def errors(self):
        self.write(self.policy, 'ci/image-stages.yml', yaml.safe_dump(self.config))
        self.write(self.policy, 'ci/release-pins.yml', yaml.safe_dump(self.pins))
        return validate_routing(self.root, self.policy)

    def test_new_task_inherits_stage_without_task_registration(self):
        self.workflow()
        self.assertEqual(self.errors(), [])

    def test_swapped_same_repository_same_default_is_rejected(self):
        self.workflow(binding='fit_docker_image')
        self.assertTrue(any('stage' in e.lower() and 'NeverRegistered' in e for e in self.errors()))

    def test_nested_workflows_scatter_and_conditional(self):
        self.workflow('workflows/child.wdl', wrapper='if (true)')
        self.workflow(child=True, wrapper='scatter (i in [1, 2])')
        self.assertEqual(self.errors(), [])

    def test_wrong_intermediate_binding(self):
        self.workflow('workflows/child.wdl')
        self.workflow(child=True, binding='fit_docker_image')
        self.assertTrue(any('stage' in e.lower() and 'child' in e for e in self.errors()))

    def test_omitted_forwarding_does_not_use_default(self):
        self.workflow(omit=True, default=True)
        self.assertTrue(any('docker_image' in e and 'forward' in e for e in self.errors()))

    def test_hard_coded_runtime_is_rejected(self):
        self.workflow(runtime='"same-image"')
        self.assertTrue(any('runtime' in e for e in self.errors()))

    def test_unknown_and_missing_scripts(self):
        for path in ('/opt/prepare_qtl/scripts/unregistered/a.R',
                     '/opt/prepare_qtl/scripts/cell_type_specific_expression/export/missing.R'):
            with self.subTest(path=path):
                self.workflow(command='Rscript ' + path)
                self.assertTrue(any(path in e for e in self.errors()))

    def test_dynamic_targets_and_inline_programs_fail_closed(self):
        for target in ('"~{script}"', '"${script}"', '"$(echo foo)"',
                       '-e "print(1)"', '-', '-m module', '$TARGET'):
            with self.subTest(target=target):
                self.workflow(command='python3 ' + target)
                self.assertTrue(any('contract' in e.lower() for e in self.errors()))

    def test_dynamic_target_error_names_the_unresolved_input(self):
        self.workflow(command='Rscript "~{script}"')
        self.assertTrue(any('~{script}' in e and 'NeverRegistered' in e for e in self.errors()))

    def test_valid_flags_quotes_continuations_and_multiple_same_stage_scripts(self):
        self.script('export/with space.py')
        self.workflow(command='Rscript --vanilla \\\n "'+ self.command.split()[1] +'" --input "~{script}"\n'
                      'python3 -u -W ignore "/opt/prepare_qtl/scripts/cell_type_specific_expression/export/with space.py"')
        self.assertEqual(self.errors(), [])
        task = WDL.load(str(self.root / 'workflows/main.wdl')).tasks[0]
        self.assertEqual(task_stages(task, self.config, self.root), {'cell_export'})

    def test_two_owned_stages_require_review(self):
        self.workflow(command=self.command + '\nRscript /opt/prepare_qtl/scripts/cell_type_specific_expression/fit/fit.R')
        self.assertTrue(any('stage' in e.lower() for e in self.errors()))

    def test_shared_helper_accepts_propagated_consumer_stage(self):
        self.workflow('workflows/child.wdl', command='Rscript /opt/prepare_qtl/scripts/common/helper.R')
        self.workflow(child=True, binding='fit_docker_image')
        self.assertEqual(self.errors(), [])

    def test_candidate_policy_cannot_reclassify_scripts(self):
        self.workflow(binding='fit_docker_image')
        self.write(self.root, 'ci/image-stages.yml', yaml.safe_dump({'stages': {}}))
        self.write(self.root, 'ci/release-pins.yml', yaml.safe_dump({'stages': {}}))
        self.assertTrue(any('stage' in e.lower() for e in self.errors()))

    def test_network_and_outside_imports_fail_without_loading(self):
        for uri in ('https://example.invalid/task.wdl', '../../outside.wdl'):
            self.write(self.root, 'workflows/main.wdl', 'version 1.0\nimport "' + uri + '" as other\nworkflow Main {}')
            self.assertTrue(any('import' in e.lower() or 'source' in e.lower() for e in self.errors()))

    def test_image_expression_is_not_evaluated(self):
        self.workflow(binding='read_string(write_lines(["bad"]))')
        self.assertTrue(any('input' in e.lower() for e in self.errors()))

    def test_unregistered_runtime_and_interpreter_wrapper_fail(self):
        for command in ('echo hello', 'bash -c "' + self.command + '"'):
            self.workflow(command=command)
            self.assertTrue(any('contract' in e.lower() for e in self.errors()))

    def test_script_words_in_echo_do_not_register_a_runtime(self):
        self.workflow(command='echo ' + self.command)
        self.assertTrue(any('contract' in e.lower() for e in self.errors()))

    def test_valid_script_does_not_hide_an_indirect_invocation(self):
        for extra in ('python3 "~{script}"', '"${runner}" script.R',
                      'bash -c "python3 ${target}"', 'eval "${command}"',
                      'MODE=test python3 "${target}"',
                      'command python3 "${target}"', 'timeout 3 python3 "${target}"'):
            with self.subTest(extra=extra):
                self.workflow(command=self.command + '\n' + extra)
                self.assertTrue(any('contract' in e.lower() for e in self.errors()))

    def test_literal_cat_heredoc_is_data(self):
        self.workflow(command="cat > paths.txt <<'PATHS'\n~{script}\nPATHS\n" + self.command)
        self.assertEqual(self.errors(), [])

    def test_valid_script_does_not_hide_timeout_shell_with_dynamic_target(self):
        self.workflow(command=self.command + '\n' +
                      "timeout 10 bash -c 'python3 \"$SCRIPT\"'")
        task = WDL.load(str(self.root / 'workflows/main.wdl')).tasks[0]
        with self.assertRaisesRegex(
                ValueError, r"workflows/main.wdl task NeverRegistered: .*bash.*contract"):
            task_stages(task, self.config, self.root)

    def test_valid_script_does_not_hide_timeout_shell_with_other_stage(self):
        self.script('fit_tca.R')
        self.config['stages']['cell_fit']['sources'] = [
            'scripts/cell_type_specific_expression/fit_tca.R']
        self.workflow(command=self.command + '\n' +
                      "timeout 10 bash -c 'Rscript /opt/prepare_qtl/scripts/"
                      "cell_type_specific_expression/fit_tca.R'")
        task = WDL.load(str(self.root / 'workflows/main.wdl')).tasks[0]
        with self.assertRaisesRegex(
                ValueError, r"workflows/main.wdl task NeverRegistered: .*bash.*contract"):
            task_stages(task, self.config, self.root)

    def test_fixed_xargs_transfer_programs_remain_supported(self):
        for count, program in enumerate(self.TRANSFER_PROGRAMS, start=2):
            with self.subTest(count=count):
                self.workflow(command=self.command + '\n' +
                              f"xargs -0 -n {count} -P 2 bash -c '{program}' _")
                self.assertEqual(self.errors(), [])

    def test_comment_quotes_do_not_change_literal_transfer_program(self):
        self.workflow(command=self.command + '\n' +
                      "# Don't change transfer settings\n" +
                      '''xargs -0 -n 2 -P 2 bash -c 'gsutil cp "$1" "$2"' _''')
        self.assertEqual(self.errors(), [])

    def test_modified_xargs_transfer_program_requires_review(self):
        for count, program in enumerate(self.TRANSFER_PROGRAMS, start=2):
            for extra in ('; python3 "$SCRIPT"', '; touch unreviewed'):
                with self.subTest(count=count, extra=extra):
                    self.workflow(command=self.command + '\n' +
                                  f"xargs -0 -n {count} -P 2 bash -c '" +
                                  program.rstrip() + extra + "' _")
                    errors = self.errors()
                    self.assertTrue(any('bash' in error and 'contract' in error
                                        for error in errors), errors)

    def test_double_quoted_transfer_program_requires_review(self):
        for count, program in enumerate(self.TRANSFER_PROGRAMS, start=2):
            with self.subTest(count=count):
                expanded = '"' + program.replace('"', '\\"') + '"'
                self.workflow(command=self.command + '\n' +
                              f'xargs -0 -n {count} -P 2 bash -c {expanded} _')
                errors = self.errors()
                self.assertTrue(any('bash' in error and 'contract' in error
                                    for error in errors), errors)

    def test_mixed_quoted_transfer_program_requires_review(self):
        self.workflow(command=self.command + '\n' +
                      '''xargs -0 -n 2 -P 2 bash -c 'gsutil cp '"\\"$1\\" \\"$2\\"" _''')
        errors = self.errors()
        self.assertTrue(any('bash' in error and 'contract' in error
                            for error in errors), errors)

    def test_ordinary_read_variable_named_source_remains_supported(self):
        self.workflow(command='read -r source destination\n' + self.command)
        self.assertEqual(self.errors(), [])

    def test_xargs_replacement_cannot_change_a_fixed_transfer_program(self):
        self.workflow(command=self.command + '\n' +
                      "xargs -I gsutil bash -c 'gsutil cp \"$1\" \"$2\"' _")
        errors = self.errors()
        self.assertTrue(any('bash' in error and 'contract' in error
                            for error in errors), errors)

    def test_nested_interpreter_substitutions_require_review(self):
        for nested in ('python3 "~{script}"',
                       'Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/fit/fit.R'):
            for substitution in ('$(' + nested + ')', '`' + nested + '`'):
                with self.subTest(substitution=substitution):
                    self.workflow(command=self.command + ' --arg "' + substitution + '"')
                    self.assertTrue(any('substitution' in e and 'contract' in e
                                        and 'NeverRegistered' in e for e in self.errors()))

    def test_current_log_substitutions_remain_supported(self):
        self.workflow(command='printf "%s\\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"\n'
                      'count="$(wc -l < outputs/data.tsv)"\n' + self.command)
        self.assertEqual(self.errors(), [])

    def test_source_symlink_cannot_escape_candidate_root(self):
        script = self.root / 'scripts/cell_type_specific_expression/export/new.R'
        script.unlink()
        outside = Path(self.temp.name) / 'outside.R'
        outside.write_text('# outside candidate')
        script.symlink_to(outside)
        self.workflow()
        self.assertTrue(any('candidate source' in e for e in self.errors()))

    def test_literal_image_binding_cannot_bypass_stage_identity(self):
        self.workflow(binding='"same-image"')
        self.assertTrue(any('unresolved stage input' in e for e in self.errors()))


if __name__ == '__main__':
    unittest.main()
