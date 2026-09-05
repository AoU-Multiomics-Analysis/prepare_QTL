#!/usr/bin/env python3
"""Resolve repository scripts to release stages without inspecting their contents."""

from fnmatch import fnmatchcase
from pathlib import Path, PurePosixPath
import subprocess


CANONICAL_RUNTIME_PREFIX = '/opt/prepare_qtl/scripts/'


def _relative_posix_path(path: str, description: str) -> str:
    if not isinstance(path, str) or not path:
        raise ValueError(f'{description} must be a nonempty string')
    if '\\' in path or path.startswith('/'):
        raise ValueError(f'{description} must be a repository-relative POSIX path: {path!r}')
    parts = path.split('/')
    if any(part in ('', '.', '..') for part in parts):
        raise ValueError(f'{description} must be normalized and cannot contain empty, . or .. components: {path!r}')
    if PurePosixPath(path).as_posix() != path:
        raise ValueError(f'{description} must be a normalized POSIX path: {path!r}')
    return path


def _runtime_path(path: str) -> str:
    if not isinstance(path, str) or not path:
        raise ValueError('Runtime script path must be a nonempty string')
    if len(path) >= 2 and path[0] == path[-1] and path[0] in ('"', "'"):
        path = path[1:-1]
    if '\\' in path or not path.startswith('/'):
        raise ValueError(f'Runtime script path must be an absolute POSIX path: {path!r}')
    parts = path.split('/')[1:]
    if any(part in ('', '.', '..') for part in parts):
        raise ValueError(f'Runtime script path must be normalized: {path!r}')
    return path


def _matches(path: str, patterns: list[str]) -> bool:
    return any(fnmatchcase(path, pattern) for pattern in patterns)


def source_patterns(stage: dict) -> list[str]:
    """Return legacy source patterns plus recursive globs for owned roots."""
    return list(stage.get('sources', [])) + [
        f'{script_root}/**' for script_root in stage.get('script_roots', [])
    ]


def stages_for_script(config: dict, source_path: str) -> set[str]:
    """Return all release stages that consume a repository-relative script."""
    source_path = _relative_posix_path(source_path, 'Source path')
    stages = {
        stage_name
        for stage_name, stage in config.get('stages', {}).items()
        if _matches(source_path, source_patterns(stage))
    }
    for shared in config.get('shared', []):
        if _matches(source_path, source_patterns(shared)):
            stages.update(shared.get('stages', []))
    return stages


def _root_covers_source(script_root: str, source_path: str) -> bool:
    return PurePosixPath(source_path).parent.as_posix() == script_root


def validate_script_roots(config: dict) -> list[str]:
    """Validate stage ownership, build coverage, and runtime path policy."""
    errors = []
    stages = config.get('stages', {})
    images = config.get('images', {})
    owned_roots = []

    for stage_name, stage in stages.items():
        for script_root in stage.get('script_roots', []):
            try:
                script_root = _relative_posix_path(script_root, 'Script root')
            except ValueError as error:
                errors.append(f'Stage {stage_name}: {error}')
                continue
            owned_roots.append((script_root, stage_name))
            image_name = stage.get('image')
            if image_name not in images:
                continue
            probe = f'{script_root}/__script_stage_probe__'
            if not _matches(probe, images[image_name].get('build_paths', [])):
                errors.append(
                    f'Script root {script_root} for stage {stage_name} is outside '
                    f'image {image_name} build paths')

    for index, (left_root, left_stage) in enumerate(owned_roots):
        for right_root, right_stage in owned_roots[index + 1:]:
            if (left_root == right_root
                    or left_root.startswith(right_root + '/')
                    or right_root.startswith(left_root + '/')):
                errors.append(
                    f'Owned script roots overlap: {left_root} ({left_stage}) and '
                    f'{right_root} ({right_stage})')

    for shared in config.get('shared', []):
        for stage_name in shared.get('stages', []):
            if stage_name not in stages:
                errors.append(f'Unknown shared script consumer: {stage_name}')

    legacy_roots = config.get('legacy_runtime_roots', [])
    for legacy_root in legacy_roots:
        try:
            legacy_root = _relative_posix_path(legacy_root, 'Legacy runtime root')
        except ValueError as error:
            errors.append(str(error))
            continue
        standard = images.get('standard')
        if standard is None:
            errors.append('Legacy runtime roots require the standard image')
        elif not _matches(f'{legacy_root}/__legacy_runtime_probe__.R',
                          standard.get('build_paths', [])):
            errors.append(
                f'Legacy runtime root {legacy_root} is outside standard image build paths')

    for runtime, source in config.get('runtime_aliases', {}).items():
        try:
            runtime = _runtime_path(runtime)
            source = _relative_posix_path(source, 'Runtime alias target')
        except ValueError as error:
            errors.append(str(error))
            continue
        if not runtime.startswith('/tmp/'):
            errors.append(f'Runtime alias must use the /tmp/ legacy prefix: {runtime}')
        if not any(_root_covers_source(root, source) for root in legacy_roots):
            errors.append(
                f'Runtime alias target is outside legacy runtime roots: {source}')
        try:
            consumers = stages_for_script(config, source)
        except ValueError as error:
            errors.append(str(error))
        else:
            if not consumers:
                errors.append(f'Runtime alias target has no registered stage: {source}')

    return errors


def _tracked_sources(root: Path) -> set[str]:
    try:
        output = subprocess.check_output(
            ['git', 'ls-files', '-z'], cwd=root, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f'Cannot list tracked sources in {root}') from error
    return {path for path in output.decode().split('\0') if path}


def source_for_runtime(config: dict, runtime_path: str, root: Path) -> str:
    """Resolve one supported runtime path to one repository-relative source."""
    runtime_path = _runtime_path(runtime_path)
    aliases = config.get('runtime_aliases', {})
    if runtime_path in aliases:
        source = _relative_posix_path(
            aliases[runtime_path], 'Runtime alias target')
        if source not in _tracked_sources(root):
            raise ValueError(f'Runtime alias target is not tracked: {source}')
        if not any(_root_covers_source(script_root, source)
                   for script_root in config.get('legacy_runtime_roots', [])):
            raise ValueError(f'Runtime alias target is outside legacy runtime roots: {source}')
        if not stages_for_script(config, source):
            raise ValueError(f'Runtime alias target has no registered stage: {source}')
        return source

    if runtime_path.startswith(CANONICAL_RUNTIME_PREFIX):
        suffix = runtime_path[len(CANONICAL_RUNTIME_PREFIX):]
        source = _relative_posix_path(f'scripts/{suffix}', 'Canonical runtime source')
        if not stages_for_script(config, source):
            raise ValueError(f'Unknown runtime script: {runtime_path}')
        return source

    if runtime_path.startswith('/tmp/'):
        basename = runtime_path.removeprefix('/tmp/')
        if '/' in basename or not basename.endswith('.R'):
            raise ValueError(f'Unknown runtime script: {runtime_path}')
        tracked = _tracked_sources(root)
        candidates = sorted(
            source for source in tracked
            if PurePosixPath(source).name == basename
            and any(_root_covers_source(script_root, source)
                    for script_root in config.get('legacy_runtime_roots', []))
            and stages_for_script(config, source)
        )
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise ValueError(
                f'Ambiguous runtime script {runtime_path}: {", ".join(candidates)}')

    raise ValueError(f'Unknown runtime script: {runtime_path}')
