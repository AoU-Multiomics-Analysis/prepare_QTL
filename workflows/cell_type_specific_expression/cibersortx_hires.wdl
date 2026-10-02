version 1.0

workflow CIBERSORTxHiRes {
  input {
    File mixture
    File signature
    File fractions
    File gene_subset
    String username
    String token
    Int threads = 8
    Int nsampling = 1
    Int nsampling2 = 1
    Int memory_gb = 16
    Int disk_gb = 20
  }

  call RunHiRes {
    input:
      mixture = mixture,
      signature = signature,
      fractions = fractions,
      gene_subset = gene_subset,
      username = username,
      token = token,
      threads = threads,
      nsampling = nsampling,
      nsampling2 = nsampling2,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    Array[File] expression_matrices = RunHiRes.expression_matrices
    File input_validation = RunHiRes.input_validation
    File output_validation = RunHiRes.output_validation
    File run_log = RunHiRes.run_log
    File task_stdout = RunHiRes.task_stdout
    File task_stderr = RunHiRes.task_stderr
  }
}

task RunHiRes {
  input {
    File mixture
    File signature
    File fractions
    File gene_subset
    String username
    String token
    Int threads = 8
    Int nsampling = 1
    Int nsampling2 = 1
    Int memory_gb = 16
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail
    {
    set -euo pipefail
    task_root="$PWD"
    log() { printf '[%s] %s\n' "$(date -u +%FT%TZ)" "$*"; }
    log "stage=CIBERSORTxHiRes start_time=$(date -u +%FT%TZ) Check localized inputs."

    # Escape apostrophes for Bash single quotes only at command rendering.
    # File inputs remain File values in the workflow and task interfaces.
    python3 - \
      --mixture '~{sub(mixture, "'", "'\"'\"'")}' \
      --signature '~{sub(signature, "'", "'\"'\"'")}' \
      --fractions '~{sub(fractions, "'", "'\"'\"'")}' \
      --gene-subset '~{sub(gene_subset, "'", "'\"'\"'")}' \
      --threads ~{threads} --nsampling ~{nsampling} --nsampling2 ~{nsampling2} \
      --memory-gb ~{memory_gb} --disk-gb ~{disk_gb} <<'PY'
    import argparse, csv, hashlib, json, math, os
    from pathlib import Path

    parser = argparse.ArgumentParser()
    for flag in ('mixture', 'signature', 'fractions', 'gene-subset'):
        parser.add_argument('--' + flag, required=True)
    for flag in ('threads', 'nsampling', 'nsampling2', 'memory-gb', 'disk-gb'):
        parser.add_argument('--' + flag, type=int, required=True)
    args = parser.parse_args()
    def require(condition, message):
        if not condition:
            raise SystemExit('Input error: ' + message)
    for flag in ('threads', 'nsampling', 'nsampling2', 'memory_gb', 'disk_gb'):
        require(getattr(args, flag) > 0, flag + ' must be positive.')
    paths = {name: Path(getattr(args, name)) for name in ('mixture', 'signature', 'fractions', 'gene_subset')}
    for name, path in paths.items():
        if '://' in str(path) or not path.is_file() or not os.access(str(path), os.R_OK):
            raise SystemExit('Localization error: ' + name + ' must be a readable local File input: ' + str(path))

    def matrix(path):
        raw = path.read_bytes()
        require(not raw.startswith(b'\xef\xbb\xbf') and b'\x00' not in raw and b'"' not in raw,
                str(path) + ' has unsupported text characters.')
        rows = list(csv.reader(raw.decode('utf-8').splitlines(), delimiter='\t'))
        require(len(rows) > 1 and len(rows[0]) > 1, str(path) + ' is empty or has no numeric columns.')
        header, ids = rows[0][1:], [row[0] if row else '' for row in rows[1:]]
        require(all(len(row) == len(rows[0]) for row in rows), str(path) + ' is not rectangular.')
        require(len(set(header)) == len(header) and len(set(ids)) == len(ids), str(path) + ' has duplicate labels.')
        require(all(label and label == label.strip() for label in header + ids), str(path) + ' has empty or padded labels.')
        try:
            values = [[float(value) for value in row[1:]] for row in rows[1:]]
        except ValueError:
            raise SystemExit('Input error: ' + str(path) + ' contains a nonnumeric value.')
        require(all(math.isfinite(value) for row in values for value in row), str(path) + ' contains a missing or nonfinite value.')
        return header, ids, values

    samples, genes, mixture_values = matrix(paths['mixture'])
    cell_types, signature_genes, signature_values = matrix(paths['signature'])
    fraction_types, fraction_samples, fraction_values = matrix(paths['fractions'])
    require(fraction_types == cell_types,
            'Use fraction-only results. Their columns must match the signature cell types in the same order. Remove P-value, Correlation and RMSE.')
    fraction_input_sample_count = len(fraction_samples)
    fraction_index = {sample: index for index, sample in enumerate(fraction_samples)}
    missing_samples = [sample for sample in samples if sample not in fraction_index]
    require(not missing_samples, 'Mixture samples missing from fractions: ' + ', '.join(missing_samples))
    # Select by sample ID, then put fraction rows in mixture column order.
    fraction_indices = [fraction_index[sample] for sample in samples]
    fraction_values = [fraction_values[index] for index in fraction_indices]
    require(len(samples) >= 4 * len(cell_types), 'HiRes needs at least four samples per cell type for the default window.')
    require(all(0 <= value <= 1 for row in fraction_values for value in row), 'Fractions must be between zero and one.')
    require(all(abs(sum(row) - 1) <= 1e-6 for row in fraction_values), 'Each fraction row must sum to one.')
    require(all(sum(row[column] for row in fraction_values) > 0 for column in range(len(cell_types))), 'Each cell type must have a positive total fraction.')
    require(set(signature_genes).issubset(genes), 'Some signature genes are absent from the mixture. Supply the shared-gene signature.')
    panel = paths['gene_subset'].read_text().splitlines()
    require(panel and len(set(panel)) == len(panel) and all(gene and gene == gene.strip() and '\t' not in gene for gene in panel),
            'The gene subset must contain one unique gene per line, with no header.')
    require(set(panel).issubset(genes), 'Some subset genes are absent from the mixture.')
    require(all(sum(row[column] for row in mixture_values) > 0 for column in range(len(samples))), 'Each mixture sample must have a positive expression total.')
    report = {
        'sample_count': len(samples), 'samples': samples, 'cell_types': cell_types,
        'fraction_input_sample_count': fraction_input_sample_count,
        'fraction_retained_sample_count': len(samples),
        'fraction_dropped_sample_count': fraction_input_sample_count - len(samples),
        'mixture_gene_count': len(genes), 'signature_gene_count': len(signature_genes), 'subset_genes': panel,
        'threads': args.threads, 'nsampling': args.nsampling, 'nsampling2': args.nsampling2,
        'negative_mixture_values': sum(value < 0 for row in mixture_values for value in row),
        'negative_signature_values': sum(value < 0 for row in signature_values for value in row),
        'sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}
    }
    Path('input_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    staging = Path('inputs').resolve()
    results = Path('results').resolve()
    staging.mkdir()
    results.mkdir()
    for name in ('mixture', 'signature', 'gene_subset'):
        (staging / (name + '.txt')).symlink_to(paths[name].resolve())
    # HiRes reads --cibresults from its output directory.
    # Preserve the header and numeric text from the original fraction rows.
    fraction_lines = paths['fractions'].read_text().splitlines()
    selected_lines = [fraction_lines[0]] + [fraction_lines[index + 1] for index in fraction_indices]
    (results / 'fractions.txt').write_text('\n'.join(selected_lines) + '\n')
    print('stage=CIBERSORTxHiRes fraction_samples=input:%d,retained:%d,dropped:%d Fractions matched to mixture samples.' %
          (fraction_input_sample_count, len(samples), fraction_input_sample_count - len(samples)), flush=True)
    for destination, source in ((Path('/src/data'), staging), (Path('/src/outdir'), results)):
        require(not destination.exists() and not destination.is_symlink(), str(destination) + ' already exists in the image.')
        destination.symlink_to(source, target_is_directory=True)
    print('stage=CIBERSORTxHiRes dimensions=samples:%d,cell_types:%d,subset_genes:%d Inputs ready.' % (len(samples), len(cell_types), len(panel)), flush=True)
    if report['negative_mixture_values'] or report['negative_signature_values']:
        print('Adjusted inputs contain negative values. Values are preserved; see input_validation.json.', flush=True)
    PY

    log 'Start HiRes: ~{threads} threads; nsampling=~{nsampling}; nsampling2=~{nsampling2}.'
    cd "/src"
    test -x ./CIBERSORTxHiRes
    # A status line confirms that the command is active. It is not a progress percentage.
    python3 -u -c 'import time; start=time.monotonic();
    while True:
        time.sleep(60)
        print("[HiRes status] Command active for %d minutes." % ((time.monotonic()-start)//60), flush=True)' &
    status_pid=$!
    trap 'kill "$status_pid" 2>/dev/null || true' EXIT
    set +e
    ./CIBERSORTxHiRes \
      --username '~{sub(username, "'", "'\"'\"'")}' \
      --token '~{sub(token, "'", "'\"'\"'")}' \
      --mixture mixture.txt --sigmatrix signature.txt --subsetgenes gene_subset.txt \
      --cibresults fractions.txt --label ZNF804A_controls \
      --QN FALSE --variableonly FALSE \
      --threads ~{threads} --nsampling ~{nsampling} --nsampling2 ~{nsampling2} \
      --outdir "/src/outdir" 2>&1 \
      | python3 -u -c 'import sys;
    for line in sys.stdin:
        if not line.startswith((">[Options] token:", ">[Options] username:")):
            print(line, end="", flush=True)'
    native_status=$?
    set -e
    kill "$status_pid" 2>/dev/null || true
    wait "$status_pid" 2>/dev/null || true
    trap - EXIT
    log "HiRes exited with status $native_status."
    if [ "$native_status" -ne 0 ]; then
      exit "$native_status"
    fi
    cd "$task_root"
    log 'Check sample-level expression outputs.'
    python3 - <<'PY'
    import csv, json, math
    from pathlib import Path
    expected = json.loads(Path('input_validation.json').read_text())
    files = sorted(Path('results').glob('CIBERSORTxHiRes*_Window*.txt'))
    remaining = set(expected['cell_types'])
    matrices = []
    for path in files:
        matches = [cell for cell in remaining if ('_' + cell + '_Window') in path.name]
        if len(matches) != 1:
            raise SystemExit('Output error: unexpected or duplicate cell type in ' + path.name)
        cell = matches[0]
        with path.open(newline='') as stream:
            rows = list(csv.reader(stream, delimiter='\t'))
        if not rows or rows[0][1:] != expected['samples']:
            raise SystemExit('Output error: sample columns do not match the mixture in ' + path.name)
        if any(len(row) != len(rows[0]) for row in rows):
            raise SystemExit('Output error: nonrectangular matrix ' + path.name)
        ids = [row[0] for row in rows[1:]]
        if len(ids) != len(set(ids)) or not set(expected['subset_genes']).issubset(ids):
            raise SystemExit('Output error: duplicate or missing subset genes in ' + path.name)
        values = [value for row in rows[1:] for value in row[1:]]
        missing = sum(value.strip().upper() in ('', 'NA', 'NAN') for value in values)
        observed = [float(value) for value in values if value.strip().upper() not in ('', 'NA', 'NAN')]
        if not all(math.isfinite(value) for value in observed):
            raise SystemExit('Output error: nonfinite estimates in ' + path.name)
        matrices.append({'cell_type': cell, 'file': str(path), 'gene_count': len(ids),
                         'missing_values': missing, 'values_equal_to_one': sum(value == 1 for value in observed)})
        remaining.remove(cell)
    if remaining:
        raise SystemExit('Output error: missing cell-type matrices: ' + ', '.join(sorted(remaining)))
    report = {'sample_count': expected['sample_count'], 'cell_types': expected['cell_types'], 'matrices': matrices,
              'scale': 'Native HiRes expression scale; units are unverified.',
              'value_one_note': 'HiRes documents 1 as an insufficient-estimation marker. Values are preserved.'}
    Path('output_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('stage=CIBERSORTxHiRes outputs=results/CIBERSORTxHiRes*_Window*.txt,input_validation.json,output_validation.json,hires.log Outputs ready: %d matrices with %d matching samples.' % (len(files), expected['sample_count']), flush=True)
    PY
    log "stage=CIBERSORTxHiRes completion_time=$(date -u +%FT%TZ) Task complete."
    } 2>&1 | tee hires.log
  >>>

  output {
    Array[File] expression_matrices = glob("results/CIBERSORTxHiRes*_Window*.txt")
    File input_validation = "input_validation.json"
    File output_validation = "output_validation.json"
    File run_log = "hires.log"
    File task_stdout = stdout()
    File task_stderr = stderr()
  }

  runtime {
    docker: "cibersortx/hires@sha256:e8da6850311d163e33a343d29a0d2ffc8b18c1ec4604995b8285c7bc2017c83e"
    cpu: threads
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} SSD"
    bootDiskSizeGb: 20
    preemptible: 0
    maxRetries: 0
  }
}
