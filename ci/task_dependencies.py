"""Check explicit task dependencies, including transitive cell-type R modules.

The module scan is conservative: identifiers in strings/comments may add an
extra dependency. Runtime loading uses the checked-in module manifests, not
source-code evaluation. New dependency patterns require an explicit review.
"""
from pathlib import Path
import re
import ast

CELL = 'scripts/cell_type_specific_expression/'
TOKEN = re.compile(r'[A-Za-z_.][A-Za-z0-9_.]*')
EXPORT = re.compile(r'^([A-Za-z_.][A-Za-z0-9_.]*)\s*<-', re.M)


def module_closure(root, scripts):
    exports = {}
    for path in sorted((root / CELL / 'R').glob('*.R')):
        for name in EXPORT.findall(path.read_text()):
            if name in exports:
                raise ValueError('Ambiguous R module symbol: ' + name)
            exports[name] = path
    seen, todo = set(), [root / p for p in scripts]
    while todo:
        for symbol in TOKEN.findall(todo.pop().read_text()):
            if symbol in exports and exports[symbol] not in seen:
                seen.add(exports[symbol])
                todo.append(exports[symbol])
    return {p.relative_to(root).as_posix() for p in seen}


def validate_task_dependencies(config, root):
    errors = []
    tasks = {}
    for name, item in config['stages'].items():
        task = item.get('task')
        if not task:  # Version-1 legacy policy is still readable.
            continue
        key = (task['path'], task['name'])
        if key in tasks:
            errors.append('Duplicate task dependency record: ' + str(key))
        tasks[key] = name
        sources = set(item.get('sources', []))
        for source in sources:
            if source.startswith('/') or '..' in Path(source).parts:
                errors.append('Unsafe task source: ' + source)
            elif not any(root.glob(source)):
                errors.append('Missing task dependency: ' + source)
        entries = [p for p in sources if p.startswith(CELL) and
                   Path(p).parent.name in ('estimation', 'fit', 'export', 'downstream')]
        for source in sorted(sources):
            source_path = root / source
            if source_path.suffix != '.py' or not source_path.is_file():
                continue
            for node in ast.walk(ast.parse(source_path.read_text())):
                if isinstance(node, ast.ImportFrom) and node.level:
                    if node.level != 1:
                        errors.append(name + ': unsupported relative Python import')
                        continue
                    modules = [node.module] if node.module else [a.name for a in node.names]
                    for module in modules:
                        dependency = source_path.parent.joinpath(*module.split('.')).with_suffix('.py')
                        relative = dependency.relative_to(root).as_posix()
                        if relative not in sources:
                            errors.append(name + ': unregistered Python helper: ' + relative)
        # Literal local R helper imports must be registered for each caller.
        for source in sorted(sources):
            source_path = root / source
            if source_path.suffix != '.R' or not source_path.is_file():
                continue
            for helper in re.findall(r"[\"']([A-Za-z0-9_.-]+[.]R)[\"']", source_path.read_text()):
                if helper == 'bootstrap.R':
                    continue
                candidates = list((root / 'scripts').rglob(helper))
                if len(candidates) == 1:
                    dependency = candidates[0].relative_to(root).as_posix()
                    if dependency not in sources:
                        errors.append(name + ': unregistered R helper: ' + dependency)
        if not entries:
            continue
        required = module_closure(root, entries)
        missing = required - sources
        if missing:
            errors.append(name + ': missing transitive R modules: ' + ', '.join(sorted(missing)))
        manifests = [p for p in sources if p.startswith(CELL + 'modules/')]
        if len(manifests) != 1:
            errors.append(name + ': expected one task module manifest')
            continue
        manifest = root / manifests[0]
        if not manifest.is_file():
            continue
        actual = manifest.read_text().splitlines()
        expected = sorted(Path(p).name for p in sources if p.startswith(CELL + 'R/'))
        if sorted(actual) != expected or len(set(actual)) != len(actual):
            errors.append(name + ': module manifest does not match task dependencies')
    return errors
