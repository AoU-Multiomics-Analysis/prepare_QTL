# Terra file paths: rules and failure record

## Default argument interface

Use typed WDL inputs, safely quoted named command-line arguments, and task-local
file lists for file arrays. Do not use a JSON wrapper by default to pass script
arguments. WDL handles localization and output collection; scripts read their
arguments. JSON remains suitable for structured data and output manifests.

`FilterCellTypeBeds` now uses named arguments and a BED-path list. The list is
created with `write_lines(cell_type_beds)` in the command placeholder, then copied
to a task-local file. Do not move it to workflow scope or wrap the generated File
in `sub`. File paths must not contain embedded newline characters in this interface.
The JSON example below records the earlier failure; it is not the current design.

## Incident: generated reference-filter config, 2026-09-05

`FilterCellTypeBeds` failed before reading the BED files. R tried to open:

```text
gs://<bucket>/<submission>/call-FilterCellTypeBeds/write_json_<hash>.tmp
```

The error was `cannot open the connection` / `No such file or directory`.
This was a path-handling defect, not a missing BED or invalid expression scale.

The command contained this expression:

```wdl
'~{sub(write_json(filter_config), "'", "'\"'\"'")}'
```

`write_json` returns a new File. In managed Cromwell, its path can be a cloud URI.
The outer `sub` converts that File to a String before final command-placeholder
file mapping. R then receives the cloud URI instead of the localized filename.

Keep the generated result File-typed:

```wdl
'~{write_json(filter_config)}'
```

The shell quotes remain outside the WDL expression. They do not change its type.
For this task, the file is a Cromwell-generated temporary name. Input filenames
remain safely serialized inside JSON rather than inserted into shell syntax.
If a different backend permits shell-special characters in generated execution
paths, verify a safe quoting strategy on that backend; do not reintroduce a
string conversion around the generated File.

## Three separate boundaries

1. **Workflow scope:** Terra's workflow engine need not have a local filesystem.
   Do not call file-writing functions in workflow declarations, call inputs,
   scatter/conditional expressions, or outputs.
2. **Incoming task Files:** preserve File types, including fields in structs.
   For configurations that contain paths to files the task must open, serialize
   during command rendering, after input localization. Cromwell does not rewrite
   cloud paths embedded in the contents of an already-written JSON or text file.
3. **Newly generated command Files:** a command-time writer returns a File that
   still needs final path mapping. Pass that result directly to the placeholder.
   Do not wrap it in `sub`, concatenation, or another string-producing operation.

Escaping an existing File identifier at command-time lookup and wrapping a newly
generated File are not equivalent: the existing identifier is mapped on lookup;
the newly generated File did not exist at that point.

Cloud URLs intentionally recorded in final manifests are metadata, not local
filenames. Preserve those URLs as strings. This exception does not authorize
passing them to local-only R readers.

## Why the old tests passed

The regression helper localized incoming cloud Files, but its writer returned
local filenames. Converting those filenames to Strings did not make the test
fail. Womtool syntax validation also cannot detect this runtime path defect.

`tests/cell_type_specific_expression/test_reference_filter_wdl.py` now checks
the named arguments and generated BED list. It verifies that the generated list
remains File-typed and maps it to a local file. It checks present and absent
reference inputs, quoted paths, and the real R commands with localized files.

These are local regression tests, not a complete Terra execution. This correction
must still be checked in a user-approved Terra run. Do not claim it was tested
on Terra based only on MiniWDL, Womtool, or container tests.

## Cromwell source references

- [Identifier lookup maps incoming File values](https://github.com/broadinstitute/cromwell/blob/87/wdl/transforms/new-base/src/main/scala/wdl/transforms/base/linking/expression/values/LookupEvaluators.scala)
- [write_json returns a File; sub returns a String](https://github.com/broadinstitute/cromwell/blob/87/wdl/transforms/new-base/src/main/scala/wdl/transforms/base/linking/expression/values/EngineFunctionEvaluators.scala)
- [File mapping preserves the distinction from strings](https://github.com/broadinstitute/cromwell/blob/87/wom/src/main/scala/wom/WomFileMapper.scala)
