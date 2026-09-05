"""Validate literal script stages and image forwarding without running WDL code.

Only direct input references carry stage identities. No WDL expression is
evaluated. Shell tokenization recognizes literal Rscript/python file targets;
inline programs and compiled tools need a named, trusted contract.
"""
from fnmatch import fnmatchcase
from pathlib import Path
import re
import shlex

import WDL
import yaml

from script_stages import source_for_runtime, stages_for_script


# These existing commands cannot be classified from a repository script path.
INLINE_EXTERNAL_TASK_CONTRACTS = {
    ('workflows/methylation/merge_methylation.wdl', 'ShardMethylationManifest'):
        ('methylation', 'Rscript', '-e'),
    ('workflows/methylation/ProcessMethylationSample.wdl', 'FilterMethylationSample'):
        ('methylation_rust', 'methylation-filter', None),
    ('workflows/methylation/cohort_aggregation.wdl', 'MergeMethylationChromosome'):
        ('methylation_rust', 'methylation-chromosome-merge', None),
}
_PLACEHOLDER = '__UNRESOLVED_WDL_PLACEHOLDER__'
_INTERPRETERS = {'Rscript', 'python', 'python3'}
_IDENTIFIER = re.compile(r'^[A-Za-z_][A-Za-z_0-9]*$')
# Existing xargs transfer programs contain no repository script invocation.
# Ignore line indentation only; changed programs need a reviewed contract.
_FIXED_TRANSFER_SHELL_PROGRAMS = {
    '2': 'gsutil cp "$1" "$2"',
    '3': '\n'.join((
        'sample_id="$1"',
        'source_path="$2"',
        'local_path="$3"',
        'if [[ "$source_path" == gs://* ]]; then',
        'gsutil -q cp "$source_path" "$local_path"',
        'else',
        'if [ ! -f "$source_path" ]; then',
        'echo "Input BED file for ${sample_id} is not accessible inside the task: ${source_path}" >&2',
        'exit 1',
        'fi',
        'cp "$source_path" "$local_path"',
        'fi',
    )),
}


def _without_literal_data(command):
    """Skip the repository's quoted cat heredoc subset, never shell programs."""
    lines = iter(command.splitlines(keepends=True))
    output = []
    for line in lines:
        match = re.fullmatch(r"\s*cat\s+>\s+[A-Za-z0-9_./-]+\s+<<'([A-Za-z0-9_]+)'\s*", line)
        if match:
            for data in lines:
                if data.rstrip('\r\n') == match[1]:
                    break
            else:
                raise ValueError('Unterminated literal cat heredoc; add a reviewed contract')
            output.append('\n')
        else:
            output.append(line)
    return ''.join(output)


def _relative(node, root):
    return Path(node.pos.abspath).resolve().relative_to(root.resolve()).as_posix()


def _check_substitutions(command):
    """Permit existing scalar log/read forms; reject all other substitutions."""
    literal_path = r'[A-Za-z0-9_./-]+'
    file_argument = rf'(?:{literal_path}|"\$[A-Za-z_][A-Za-z_0-9]*")'
    safe_forms = (
        r"date -u (?:\+%Y-%m-%dT%H:%M:%SZ|'\+%Y-%m-%dT%H:%M:%SZ')",
        r"date '\+%b %d %H:%M:%S'",
        rf'wc -l < {literal_path}',
        rf"awk 'END {{ print NR (?:- 1|\+ 0) }}' {file_argument}",
        rf"zgrep -m 1 '\^' {file_argument}",
        rf'<(?:{literal_path}|"{_PLACEHOLDER}[0-9]+__")',
    )
    # The current RNA-SeQC log uses this integer arithmetic expression.
    remainder = re.sub(rf'\$\(\({_PLACEHOLDER}[0-9]+__ \+ 1\)\)', '', command)
    for match in re.finditer(r'\$\(([^()]*)\)', remainder):
        if not any(re.fullmatch(form, match[1]) for form in safe_forms):
            raise ValueError(f'unsupported command substitution {match[0]!r}; add a reviewed contract')
    remainder = re.sub(r'\$\([^()]*\)', '', remainder)
    if '$(' in remainder or '`' in remainder:
        raise ValueError('unsupported nested or backtick command substitution; add a reviewed contract')


def task_stages(task, config: dict, source_root: Path) -> set[str]:
    """Return stages that can run every literal script in this task.

    Raise ValueError for unresolved commands, absent files, or incompatible
    stages. Multiple consumer sets are intersected, including shared helpers.
    """
    label = f'{_relative(task, source_root)} task {task.name}'
    contract = INLINE_EXTERNAL_TASK_CONTRACTS.get((_relative(task, source_root), task.name))
    placeholders = {f'{_PLACEHOLDER}{i}__': str(part)
                    for i, part in enumerate(task.command.parts) if not isinstance(part, str)}

    def display(token):
        for marker, expression in placeholders.items():
            token = token.replace(marker, expression)
        return token

    command = ''.join(part if isinstance(part, str) else f'{_PLACEHOLDER}{i}__'
                      for i, part in enumerate(task.command.parts)).replace('\\\n', '')
    command = _without_literal_data(command)
    try:
        _check_substitutions(command)
    except ValueError as error:
        raise ValueError(f'{label}: {display(str(error))}') from error
    lexer = shlex.shlex(command, posix=True, punctuation_chars=';&|()')
    lexer.whitespace = ' \t\r'
    lexer.wordchars += '$'
    try:
        tokens = list(lexer)
    except ValueError as error:
        raise ValueError(f'{label}: unsupported command; add a reviewed contract: {error}') from error
    consumers = []
    matched_contract = False
    command_start = True
    command_index = 0
    for index, token in enumerate(tokens):
        executable = Path(token).name
        if token in {'\n', ';', '&&', '||', '|', '(', 'then', 'do', 'else', 'if', 'elif'}:
            command_start = True
            continue
        if not command_start:
            if executable in {'bash', 'sh'}:
                program = tokens[index + 2] if index + 2 < len(tokens) else ''
                program = '\n'.join(line.strip() for line in program.strip().splitlines())
                prefix = tokens[command_index:index]
                if (len(prefix) == 6 and prefix[:3] == ['xargs', '-0', '-n']
                        and prefix[3] in _FIXED_TRANSFER_SHELL_PROGRAMS
                        and prefix[4] == '-P'
                        and re.fullmatch(rf'(?:[0-9]+|{_PLACEHOLDER}[0-9]+__)', prefix[5])
                        and token == 'bash'
                        and tokens[index + 1:index + 2] == ['-c']
                        and program == _FIXED_TRANSFER_SHELL_PROGRAMS[prefix[3]]):
                    continue
                raise ValueError(f'{label}: unsupported indirect command {display(token)!r}; '
                                 'add a reviewed contract')
            if executable in _INTERPRETERS:
                raise ValueError(f'{label}: interpreter {token!r} is not a direct command; '
                                 'add a reviewed contract')
            continue
        if re.match(r'^[A-Za-z_][A-Za-z_0-9]*=', token):
            continue
        command_start = False
        command_index = index
        if executable in {'bash', 'sh', 'eval', 'source', '.', 'exec', 'env'} or (
                any(mark in token for mark in (_PLACEHOLDER, '$', '`'))
                and '=' not in token):
            raise ValueError(f'{label}: unsupported indirect command {display(token)!r}; '
                             'add a reviewed contract')
        if contract and token == contract[1] and contract[2] is None:
            matched_contract = True
            consumers.append({contract[0]})
        if executable not in _INTERPRETERS:
            continue
        cursor = index + 1
        while cursor < len(tokens):
            option = tokens[cursor]
            if option == '--':
                cursor += 1
                break
            if executable == 'Rscript' and (option in {
                    '--vanilla', '--no-save', '--no-restore', '--no-environ',
                    '--no-site-file', '--no-init-file', '--verbose'}
                    or option.startswith('--default-packages=')):
                cursor += 1
            elif executable in {'python', 'python3'} and option in {
                    '-u', '-B', '-E', '-I', '-s', '-S', '-O', '-OO', '-q'}:
                cursor += 1
            elif executable in {'python', 'python3'} and option in {'-W', '-X'}:
                cursor += 2
            else:
                break
        target = tokens[cursor] if cursor < len(tokens) else '<missing>'
        if contract and token == contract[1] and target == contract[2]:
            matched_contract = True
            consumers.append({contract[0]})
            continue
        if (not target.startswith('/') or any(mark in target for mark in
                (_PLACEHOLDER, '$', '`', '*', '?', '[', '\n'))):
            raise ValueError(f'{label}: unsupported dynamic or inline interpreter target '
                             f'{display(target)!r}; add a reviewed contract')
        try:
            source = source_for_runtime(config, target, source_root)
            resolved = (source_root / source).resolve()
            if not resolved.is_relative_to(source_root.resolve()) or not resolved.is_file():
                raise ValueError(f'Runtime script {target} has no file inside candidate source: {source}')
            consumers.append(stages_for_script(config, source))
        except ValueError as error:
            raise ValueError(f'{label}: {error}') from error
    if contract and not matched_contract:
        raise ValueError(f'{label}: command no longer matches its reviewed contract')
    if not consumers:
        raise ValueError(f'{label}: unregistered runtime command; add a reviewed contract')
    allowed = set.intersection(*consumers)
    if not allowed:
        raise ValueError(f'{label}: scripts require distinct stages {consumers}; review task stages')
    return allowed


def _calls(nodes):
    for node in nodes:
        if isinstance(node, WDL.Tree.Call):
            yield node
        elif isinstance(node, (WDL.Tree.Scatter, WDL.Tree.Conditional)):
            yield from _calls(node.body)


def _runtime_input(task):
    expression = task.runtime.get('docker')
    name = str(expression)
    declarations = {d.name: d for d in task.inputs or []}
    if (not _IDENTIFIER.fullmatch(name) or name not in declarations
            or str(declarations[name].type) != 'String'):
        raise ValueError(f'task {task.name}: runtime docker {expression} must use a declared String input')
    return name


def validate_routing(source_root: Path, policy_root: Path) -> list[str]:
    """Check candidate source against policy loaded only from policy_root."""
    source_root, policy_root = source_root.resolve(), policy_root.resolve()
    config = yaml.safe_load((policy_root / 'ci/image-stages.yml').read_text())
    pins = yaml.safe_load((policy_root / 'ci/release-pins.yml').read_text())
    identities = {}
    errors = []
    for stage, targets in pins['stages'].items():
        for target in targets:
            key = (target['path'], target.get('scope', 'workflow'),
                   target.get('task'), target['input'])
            if key in identities and identities[key] != stage:
                errors.append(f'{key}: conflicting registered stage inputs')
            identities[key] = stage

    def external(path):
        return any(group.get('policy') == 'external_images_manual'
                   and any(fnmatchcase(path, pattern) for pattern in group['paths'])
                   for group in config.get('workflow_groups', []))

    async def read_local(uri, search_path, importer):
        if '://' in uri:
            raise ValueError(f'Network import is unsupported: {uri}')
        base = Path(importer.pos.abspath).parent if importer else source_root
        path = (base / uri).resolve()
        if not path.is_relative_to(source_root):
            raise ValueError(f'Import outside candidate source is unsupported: {uri}')
        return WDL.Tree.ReadSourceResult(path.read_text(), str(path))

    documents = []
    for path in sorted((source_root / 'workflows').rglob('*.wdl')):
        if external(path.relative_to(source_root).as_posix()):
            continue
        try:
            documents.append(WDL.load(str(path), read_source=read_local))
        except Exception as error:
            errors.append(f'{path.relative_to(source_root)}: cannot load source/import: {error}')

    task_info = {}

    def inspect_task(task):
        key = (_relative(task, source_root), task.name)
        if key in task_info:
            return task_info[key]
        if external(key[0]):
            task_info[key] = None
            return None
        try:
            allowed = task_stages(task, config, source_root)
            image_input = _runtime_input(task)
            declared_stage = identities.get((key[0], 'task', task.name, image_input))
            if declared_stage and declared_stage not in allowed:
                raise ValueError(f'task {task.name}: default input {image_input} stage '
                                 f'{declared_stage} does not match script stages {sorted(allowed)}')
            task_info[key] = (image_input, allowed)
        except ValueError as error:
            errors.append(f'{key[0]}: {error}')
            task_info[key] = None
        return task_info[key]

    def workflow_inputs(workflow):
        path = _relative(workflow, source_root)
        return {d.name: identities[(path, 'workflow', None, d.name)]
                for d in workflow.inputs or []
                if (path, 'workflow', None, d.name) in identities}

    def forwarded(call, image_input, environment, location):
        expression = call.inputs.get(image_input)
        if expression is None:
            errors.append(f'{location}: input {image_input} must forward a stage image')
            return None
        name = str(expression)
        if not _IDENTIFIER.fullmatch(name) or name not in environment:
            errors.append(f'{location}: unresolved stage input {image_input} = {expression}; '
                          'use a registered image input directly')
            return None
        return environment[name]

    def visit(workflow, environment, trail):
        for call in _calls(workflow.body):
            location = f'{trail} -> {_relative(call, source_root)} call {call.name}'
            if isinstance(call.callee, WDL.Tree.Task):
                info = inspect_task(call.callee)
                if info is None:
                    continue
                image_input, allowed = info
                actual = forwarded(call, image_input, environment, location)
                if actual is not None and actual not in allowed:
                    errors.append(f'{location} task {call.callee.name}: input {image_input} '
                                  f'has stage {actual}; scripts require {sorted(allowed)}')
            else:
                child_environment = {}
                for image_input in workflow_inputs(call.callee):
                    stage = forwarded(call, image_input, environment, location)
                    if stage is not None:
                        child_environment[image_input] = stage
                visit(call.callee, child_environment, location)

    for document in documents:
        for task in document.tasks:
            inspect_task(task)
        if document.workflow:
            workflow = document.workflow
            visit(workflow, workflow_inputs(workflow),
                  f'{_relative(workflow, source_root)} workflow {workflow.name}')
    return list(dict.fromkeys(errors))
