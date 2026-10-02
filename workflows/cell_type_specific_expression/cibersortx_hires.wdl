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
    Int genes_per_chunk = 1000
    Int merge_memory_gb = 4
    Int merge_disk_gb = 20
  }

  call SplitGeneSubset {
    input:
      gene_subset = gene_subset,
      genes_per_chunk = genes_per_chunk
  }

  scatter (chunk_genes in SplitGeneSubset.chunk_gene_subsets) {
    call RunHiRes {
      input:
        mixture = mixture,
        signature = signature,
        fractions = fractions,
        gene_subset = chunk_genes,
        username = username,
        token = token,
        threads = threads,
        nsampling = nsampling,
        nsampling2 = nsampling2,
        memory_gb = memory_gb,
        disk_gb = disk_gb
    }
  }

  call MergeHiRes {
    input:
      gene_subset = gene_subset,
      chunk_gene_subsets = SplitGeneSubset.chunk_gene_subsets,
      chunk_expression_matrices = flatten(RunHiRes.expression_matrices),
      input_validations = RunHiRes.input_validation,
      run_logs = RunHiRes.run_log,
      memory_gb = merge_memory_gb,
      disk_gb = merge_disk_gb
  }

  output {
    Array[File] expression_matrices = MergeHiRes.expression_matrices
    File input_validation = MergeHiRes.input_validation
    File output_validation = MergeHiRes.output_validation
    File run_log = MergeHiRes.run_log
    File task_stdout = MergeHiRes.task_stdout
    File task_stderr = MergeHiRes.task_stderr
    File chunk_plan = SplitGeneSubset.chunk_plan
    Array[File] chunk_gene_lists = SplitGeneSubset.chunk_gene_subsets
    Array[File] chunk_input_validations = RunHiRes.input_validation
    Array[File] chunk_output_validations = RunHiRes.output_validation
    Array[File] chunk_logs = RunHiRes.run_log
    Array[File] chunk_stdout = RunHiRes.task_stdout
    Array[File] chunk_stderr = RunHiRes.task_stderr
    File split_log = SplitGeneSubset.run_log
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

task SplitGeneSubset {
  input {
    File gene_subset
    Int genes_per_chunk = 1000
  }

  command <<<
    set -euo pipefail
    {
    set -euo pipefail
    printf '[%s] stage=SplitGeneSubset start_time=%s Check gene list.\n' "$(date -u +%FT%TZ)" "$(date -u +%FT%TZ)"
    python3 - --gene-subset '~{sub(gene_subset, "'", "'\"'\"'")}' --genes-per-chunk ~{genes_per_chunk} <<'PY'
    import argparse, hashlib, json, os
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-subset', required=True)
    parser.add_argument('--genes-per-chunk', type=int, required=True)
    args = parser.parse_args()
    if args.genes_per_chunk <= 0:
        raise SystemExit('Input error: genes_per_chunk must be positive.')
    path = Path(args.gene_subset)
    if '://' in args.gene_subset or not path.is_file() or not os.access(str(path), os.R_OK):
        raise SystemExit('Localization error: gene_subset must be a readable local File input: ' + args.gene_subset)
    raw = path.read_bytes()
    try:
        genes = raw.decode('utf-8').splitlines()
    except UnicodeDecodeError:
        raise SystemExit('Input error: gene_subset must be UTF-8 text.')
    if (raw.startswith(b'\xef\xbb\xbf') or b'\x00' in raw or not genes or len(set(genes)) != len(genes)
            or any(not gene or gene != gene.strip() or '\t' in gene for gene in genes)):
        raise SystemExit('Input error: gene_subset must contain one unique gene per line, with no header or blank lines.')
    Path('chunks').mkdir()
    chunks = []
    for index, start in enumerate(range(0, len(genes), args.genes_per_chunk), 1):
        chunk = genes[start:start + args.genes_per_chunk]
        chunk_path = Path('chunks/genes_%08d.txt' % index)
        chunk_path.write_text('\n'.join(chunk) + '\n', encoding='utf-8')
        chunks.append({'index': index, 'file': str(chunk_path), 'gene_count': len(chunk),
                       'sha256': hashlib.sha256(chunk_path.read_bytes()).hexdigest()})
    report = {'gene_count': len(genes), 'genes_per_chunk': args.genes_per_chunk,
              'chunk_count': len(chunks), 'gene_subset_sha256': hashlib.sha256(raw).hexdigest(), 'chunks': chunks}
    Path('chunk_plan.json').write_text(json.dumps(report, indent=2) + '\n')
    print('stage=SplitGeneSubset dimensions=genes:%d,chunks:%d,genes_per_chunk:%d' %
          (len(genes), len(chunks), args.genes_per_chunk), flush=True)
    print('stage=SplitGeneSubset outputs=chunks/genes_*.txt,chunk_plan.json,split.log Gene lists ready.', flush=True)
    PY
    printf '[%s] stage=SplitGeneSubset completion_time=%s Task complete.\n' "$(date -u +%FT%TZ)" "$(date -u +%FT%TZ)"
    } 2>&1 | tee split.log
  >>>

  output {
    Array[File] chunk_gene_subsets = glob("chunks/genes_*.txt")
    File chunk_plan = "chunk_plan.json"
    File run_log = "split.log"
  }

  runtime {
    docker: "cibersortx/hires@sha256:e8da6850311d163e33a343d29a0d2ffc8b18c1ec4604995b8285c7bc2017c83e"
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 1 SSD"
    bootDiskSizeGb: 20
    preemptible: 0
    maxRetries: 0
  }
}

task MergeHiRes {
  input {
    File gene_subset
    Array[File] chunk_gene_subsets
    Array[File] chunk_expression_matrices
    Array[File] input_validations
    Array[File] run_logs
    Int memory_gb = 4
    Int disk_gb = 20
  }

  command <<<
    set -euo pipefail
    {
    set -euo pipefail
    printf '[%s] stage=MergeHiRes start_time=%s Check localized chunks.\n' "$(date -u +%FT%TZ)" "$(date -u +%FT%TZ)"
    # Each generated list remains a File until final command-path mapping.
    # The lists contain already-localized incoming File paths, in array order.
    cat "~{write_lines(chunk_gene_subsets)}" > chunk_paths.txt
    cat "~{write_lines(chunk_expression_matrices)}" > matrix_paths.txt
    cat "~{write_lines(input_validations)}" > validation_paths.txt
    cat "~{write_lines(run_logs)}" > log_paths.txt
    python3 - \
      --gene-subset '~{sub(gene_subset, "'", "'\"'\"'")}' \
      --chunk-files chunk_paths.txt --matrix-files matrix_paths.txt \
      --validation-files validation_paths.txt --log-files log_paths.txt \
      --memory-gb ~{memory_gb} --disk-gb ~{disk_gb} <<'PY'
    import argparse, hashlib, json, math, os, re, sqlite3
    from pathlib import Path

    parser = argparse.ArgumentParser()
    for name in ('gene-subset', 'chunk-files', 'matrix-files', 'validation-files', 'log-files'):
        parser.add_argument('--' + name, required=True)
    for name in ('memory-gb', 'disk-gb'):
        parser.add_argument('--' + name, type=int, required=True)
    args = parser.parse_args()
    def require(condition, message):
        if not condition:
            raise SystemExit('Merge error: ' + message)
    require(args.memory_gb > 0 and args.disk_gb > 0, 'memory_gb and disk_gb must be positive.')
    def local_file(value):
        path = Path(value)
        if '://' in value or not path.is_file() or not os.access(str(path), os.R_OK):
            raise SystemExit('Localization error: expected a readable local File input: ' + value)
        return path
    def file_list(value):
        entries = local_file(value).read_text().splitlines()
        require(entries and all(entries), 'A required File array is empty or has a blank path.')
        return [local_file(entry) for entry in entries]
    def gene_list(path):
        raw = path.read_bytes()
        try:
            genes = raw.decode('utf-8').splitlines()
        except UnicodeDecodeError:
            raise SystemExit('Merge error: gene list must be UTF-8: ' + str(path))
        require(not raw.startswith(b'\xef\xbb\xbf') and b'\x00' not in raw and genes
                and len(genes) == len(set(genes))
                and all(gene and gene == gene.strip() and '\t' not in gene for gene in genes),
                'Invalid or duplicate gene list: ' + str(path))
        return genes
    original_path = local_file(args.gene_subset)
    genes = gene_list(original_path)
    positions = {gene: index for index, gene in enumerate(genes)}
    chunks = file_list(args.chunk_files)
    matrix_files = file_list(args.matrix_files)
    validation_files = file_list(args.validation_files)
    logs = file_list(args.log_files)
    require(len(chunks) == len(validation_files) == len(logs), 'Chunk, validation, and log counts differ.')
    chunk_genes = [gene_list(path) for path in chunks]
    all_chunk_genes = [gene for chunk in chunk_genes for gene in chunk]
    require(len(all_chunk_genes) == len(set(all_chunk_genes)), 'Genes overlap between chunks.')
    require(set(all_chunk_genes) == set(genes), 'Chunk coverage has missing or unexpected genes.')
    reports = [json.loads(path.read_text()) for path in validation_files]
    expected = reports[0]
    samples, cells = expected['samples'], expected['cell_types']
    require(samples and len(samples) == len(set(samples)) and expected['sample_count'] == len(samples), 'Invalid sample labels or count.')
    require(cells and len(cells) == len(set(cells)) and all(
        cell and not any(char in cell for char in ('/', '\\', '\r', '\n', '\x00')) for cell in cells), 'Invalid cell-type labels.')
    shared_fields = ('samples', 'sample_count', 'cell_types', 'mixture_gene_count', 'signature_gene_count',
                     'threads', 'nsampling', 'nsampling2')
    for index, report in enumerate(reports):
        require(all(report.get(key) == expected.get(key) for key in shared_fields), 'Sample, cell-type, background or sampling settings differ between chunks.')
        require(all(report['sha256'].get(key) == expected['sha256'].get(key) and report['sha256'].get(key)
                    for key in ('mixture', 'signature', 'fractions')), 'Full-input checksums differ between chunks.')
        require(report['subset_genes'] == chunk_genes[index], 'Validation genes do not match their chunk list.')
        require(report['sha256']['gene_subset'] == hashlib.sha256(chunks[index].read_bytes()).hexdigest(), 'Chunk checksum does not match validation.')
    require(len(matrix_files) == len(chunks) * len(cells), 'Expected one matrix per cell type per chunk.')
    # WDL flatten preserves chunk-major order. Basenames alone are not unique.
    by_chunk = []
    window = None
    for index in range(len(chunks)):
        group = {}
        for path in matrix_files[index * len(cells):(index + 1) * len(cells)]:
            matches = []
            for cell in cells:
                match = re.search('_' + re.escape(cell) + r'_Window([0-9]+)\.txt$', path.name)
                if match:
                    matches.append((cell, int(match.group(1))))
            # Use the complete longest label if cell names share a suffix.
            require(matches, 'Unexpected cell-type matrix: ' + path.name)
            cell, size = max(matches, key=lambda match: len(match[0]))
            require(cell not in group, 'Duplicate cell-type matrix in a chunk: ' + cell)
            require(size > 0 and (window is None or window == size), 'Window sizes differ between chunks or cell types.')
            window = size
            group[cell] = path
        require(set(group) == set(cells), 'Missing cell-type matrix in a chunk.')
        by_chunk.append(group)
    print('stage=MergeHiRes dimensions=genes:%d,samples:%d,cell_types:%d,chunks:%d' %
          (len(genes), len(samples), len(cells), len(chunks)), flush=True)
    Path('results').mkdir()
    matrices = []
    for cell in cells:
        db_path = Path('merge_rows.sqlite')
        database = sqlite3.connect(str(db_path))
        database.execute('PRAGMA cache_size=-8192')
        database.execute('CREATE TABLE rows (position INTEGER PRIMARY KEY, line TEXT NOT NULL)')
        header = None
        missing = ones = 0
        for index, group in enumerate(by_chunk):
            seen = set()
            allowed = set(chunk_genes[index])
            path = group[cell]
            with path.open(encoding='utf-8', newline='') as stream:
                current_header = stream.readline().rstrip('\r\n').split('\t')
                require(current_header[1:] == samples and (header is None or current_header == header), 'Sample columns or header differ: ' + str(path))
                header = current_header
                for line in stream:
                    text = line.rstrip('\r\n')
                    row = text.split('\t')
                    require(len(row) == len(header), 'Nonrectangular matrix: ' + str(path))
                    gene = row[0]
                    require(gene in allowed and gene not in seen, 'Unexpected or duplicate gene in matrix: ' + gene)
                    seen.add(gene)
                    for value in row[1:]:
                        if value.strip().upper() in ('', 'NA', 'NAN'):
                            missing += 1
                        else:
                            try:
                                numeric = float(value)
                            except ValueError:
                                raise SystemExit('Merge error: Nonnumeric estimate in ' + str(path))
                            require(math.isfinite(numeric), 'Nonfinite estimate in ' + str(path))
                            ones += numeric == 1
                    database.execute('INSERT INTO rows VALUES (?, ?)', (positions[gene], text))
            require(seen == allowed, 'Missing genes in matrix: ' + str(path))
            database.commit()
        output = Path('results/CIBERSORTxHiRes_merged_' + cell + '_Window%d.txt' % window)
        with output.open('w', encoding='utf-8', newline='') as stream:
            stream.write('\t'.join(header) + '\n')
            for row in database.execute('SELECT line FROM rows ORDER BY position'):
                stream.write(row[0] + '\n')
        database.close()
        db_path.unlink()
        matrices.append({'cell_type': cell, 'file': str(output), 'gene_count': len(genes),
                         'missing_values': missing, 'values_equal_to_one': ones})
        print('stage=MergeHiRes Merged cell type: ' + cell, flush=True)
    merged_input = dict(expected)
    merged_input['sha256'] = dict(expected['sha256'])
    merged_input['sha256']['gene_subset'] = hashlib.sha256(original_path.read_bytes()).hexdigest()
    merged_input.update(subset_genes=genes, chunk_count=len(chunks), chunk_gene_counts=[len(chunk) for chunk in chunk_genes])
    Path('input_validation.json').write_text(json.dumps(merged_input, indent=2) + '\n')
    output_report = {'sample_count': len(samples), 'cell_types': cells, 'matrices': matrices,
                     'chunk_count': len(chunks), 'window_size': window,
                     'scale': 'Native HiRes expression scale; units are unverified.',
                     'value_one_note': 'HiRes documents 1 as an insufficient-estimation marker. Values are preserved.'}
    Path('output_validation.json').write_text(json.dumps(output_report, indent=2) + '\n')
    # These File paths are explicit task inputs, not paths read from a manifest.
    with Path('chunk_native_logs.txt').open('w') as destination:
        for index, path in enumerate(logs, 1):
            destination.write('\n--- HiRes chunk %d ---\n' % index)
            with path.open() as source:
                for line in source:
                    destination.write(line)
    print('stage=MergeHiRes outputs=results/CIBERSORTxHiRes*_Window*.txt,input_validation.json,output_validation.json,hires.log Matrices ready.', flush=True)
    PY
    cat chunk_native_logs.txt
    printf '[%s] stage=MergeHiRes completion_time=%s Task complete.\n' "$(date -u +%FT%TZ)" "$(date -u +%FT%TZ)"
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
    cpu: 1
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} SSD"
    bootDiskSizeGb: 20
    preemptible: 0
    maxRetries: 0
  }
}
