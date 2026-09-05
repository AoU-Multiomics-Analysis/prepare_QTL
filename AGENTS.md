# Terra WDL maintenance

These pipelines target Terra's managed Cromwell environment. Use WDL 1.0.
Before changing WDL file handling, read [Terra file-path rules](docs/terra-file-paths.md).

- Do not create files at workflow scope.
- Do not use JSON configuration files as the default way to pass arguments
  from WDL tasks to scripts. Use typed WDL inputs and safely quoted named CLI
  arguments. For arrays of files, use a task-local newline-delimited file list.
  Add optional CLI arguments only when their inputs are present.
- WDL owns file localization, dependencies, and output collection. Keep JSON
  for actual structured data or output manifests; justify any exception that
  uses JSON as an argument wrapper before implementing it.
- Keep files that tasks must open File-typed until localization.
- Do not wrap a command-time `write_json`, `write_lines`, or other generated
  File result in `sub`, string concatenation, or another string conversion.
  Cromwell still needs to localize that newly generated File.
- Keep cloud paths used only as manifest metadata separate from local input paths.
- Tests must cover both incoming cloud Files and newly generated cloud Files.
- Keep command logging. Do not build Docker images locally for smoke tests;
  use GitHub Actions. Do not submit Terra jobs without user approval.
- State separately whether syntax checks, local task tests, and a complete
  Terra workflow run have passed. Local tests do not establish Terra success.
