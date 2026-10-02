version 1.0

workflow CIBERSORTxMarkers {
  input {
    File reference
    String username
    String token
    Int replicates = 5
    Float sampling = 1.0
    Float fraction = 0.0
    Int min_genes = 300
    Int max_genes = 500
    Float q_value = 0.01
    Int max_condition_number = 999
    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 100
  }

  call DeriveMarkers {
    input:
      reference = reference,
      username = username,
      token = token,
      replicates = replicates,
      sampling = sampling,
      fraction = fraction,
      min_genes = min_genes,
      max_genes = max_genes,
      q_value = q_value,
      max_condition_number = max_condition_number,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File signature_matrix = DeriveMarkers.signature_matrix
    File source_geps = DeriveMarkers.source_geps
    File reference_sample = DeriveMarkers.reference_sample
    File phenotype_classes = DeriveMarkers.phenotype_classes
    Array[File] signature_heatmaps = DeriveMarkers.signature_heatmaps
    File input_validation = DeriveMarkers.input_validation
    File output_validation = DeriveMarkers.output_validation
    File run_log = DeriveMarkers.run_log
    File task_stdout = DeriveMarkers.task_stdout
    File task_stderr = DeriveMarkers.task_stderr
  }
}

task DeriveMarkers {
  input {
    File reference
    String username
    String token
    Int replicates = 5
    Float sampling = 1.0
    Float fraction = 0.0
    Int min_genes = 300
    Int max_genes = 500
    Float q_value = 0.01
    Int max_condition_number = 999
    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 100
  }

  command <<<
    set -euo pipefail
    # Escape existing File values only during command rendering, after localization.
    # No file is created or converted to a String at workflow scope.
    python3 - \
      --reference '~{sub(reference, "'", "'\"'\"'")}' \
      --username '~{sub(username, "'", "'\"'\"'")}' \
      --token '~{sub(token, "'", "'\"'\"'")}' \
      --replicates ~{replicates} --sampling ~{sampling} --fraction ~{fraction} \
      --min-genes ~{min_genes} --max-genes ~{max_genes} --q-value ~{q_value} \
      --max-condition-number ~{max_condition_number} \
      --cpu ~{cpu} --memory-gb ~{memory_gb} --disk-gb ~{disk_gb} <<'PY'
    import argparse, collections, csv, datetime, hashlib, json, math, os, shutil
    import subprocess, sys, threading
    from pathlib import Path

    parser = argparse.ArgumentParser()
    for name in ('reference', 'username', 'token'):
        parser.add_argument('--' + name, required=True)
    for name in ('replicates', 'min-genes', 'max-genes', 'max-condition-number', 'cpu', 'memory-gb', 'disk-gb'):
        parser.add_argument('--' + name, type=int, required=True)
    for name in ('sampling', 'fraction', 'q-value'):
        parser.add_argument('--' + name, type=float, required=True)
    args = parser.parse_args()
    root = Path.cwd()
    logfile = (root / 'cibersortx_markers.log').open('w', buffering=1)
    secrets = sorted(set([args.username, args.token] + args.username.splitlines() + args.token.splitlines()), key=len, reverse=True)
    def emit(message):
        for secret in secrets:
            if secret:
                message = message.replace(secret, '[redacted]')
        print(message, flush=True)
        logfile.write(message + '\n')
    def now():
        return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    def log(message):
        emit('[%s] stage=CIBERSORTxMarkers %s' % (now(), message))
    def require(condition, message):
        if not condition:
            raise ValueError(message)
    def digest(path):
        result = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                result.update(block)
        return result.hexdigest()
    def matrix(path, context, repeated_labels=False):
        # Stream rows: a single-cell reference can contain tens of thousands of columns.
        with path.open(newline='') as stream:
            reader = csv.reader(stream, delimiter='\t')
            header = next(reader, [])
            require(len(header) > 1, context + ': matrix has no expression columns.')
            columns = header[1:]
            require(all(x and x == x.strip() for x in columns), context + ': empty or padded column label.')
            if not repeated_labels:
                require(len(set(columns)) == len(columns), context + ': duplicate column labels.')
            genes = set()
            totals = [0.0] * len(columns)
            for row in reader:
                require(len(row) == len(header), context + ': matrix is not rectangular.')
                require(row[0] and row[0] == row[0].strip() and row[0] not in genes,
                        context + ': empty, padded, or duplicate gene label.')
                genes.add(row[0])
                try:
                    values = [float(value) for value in row[1:]]
                except ValueError:
                    raise ValueError(context + ': nonnumeric expression.')
                require(all(math.isfinite(value) and value >= 0 for value in values),
                        context + ': expression must be finite and nonnegative; supply linear expression.')
                totals = [a + b for a, b in zip(totals, values)]
            require(bool(genes), context + ': matrix has no gene rows.')
            require(all(math.isfinite(value) and value > 0 for value in totals),
                    context + ': each column must have positive finite total expression.')
        return columns, genes
    def one_output(pattern):
        paths = list(results.glob(pattern))
        require(len(paths) == 1 and paths[0].is_file() and paths[0].stat().st_size > 0,
                'Output error: expected one nonempty ' + pattern)
        return paths[0]

    try:
        log('start_time=%s Check localized reference.' % now())
        require('://' not in args.reference, 'Localization error: reference is an unresolved cloud URI.')
        reference = Path(args.reference).resolve()
        require(reference.is_file() and os.access(str(reference), os.R_OK),
                'Localization error: reference is missing or unreadable.')
        require(bool(args.username.strip()) and bool(args.token.strip()), 'Input error: username and token are required.')
        require(args.replicates >= 2, 'Input error: replicates must be at least 2.')
        require(0 < args.sampling <= 1, 'Input error: sampling must be greater than 0 and at most 1.')
        require(0 <= args.fraction <= 1, 'Input error: fraction must be between 0 and 1.')
        require(0 < args.q_value <= 1, 'Input error: q_value must be greater than 0 and at most 1.')
        require(0 < args.min_genes <= args.max_genes, 'Input error: require 0 < min_genes <= max_genes.')
        require(all(getattr(args, name) > 0 for name in ('max_condition_number', 'cpu', 'memory_gb', 'disk_gb')),
                'Input error: condition number, CPU, memory, and disk requests must be positive.')
        columns, input_genes = matrix(reference, 'Input error: reference', repeated_labels=True)
        counts = dict(collections.Counter(columns))
        require(len(counts) >= 2, 'Input error: reference needs at least two cell types.')
        require(all(value >= 3 for value in counts.values()), 'Input error: each cell type needs at least three cells.')
        parameters = {name: getattr(args, name) for name in (
            'replicates', 'sampling', 'fraction', 'min_genes', 'max_genes', 'q_value', 'max_condition_number',
            'cpu', 'memory_gb', 'disk_gb')}
        input_report = {'gene_count': len(input_genes), 'cell_count': len(columns), 'cell_types': list(counts),
                        'cells_per_type': counts, 'parameters': parameters, 'sha256': digest(reference)}
        (root / 'input_validation.json').write_text(json.dumps(input_report, indent=2) + '\n')
        staging = root / 'inputs'
        results = root / 'results'
        staging.mkdir()
        results.mkdir()
        (staging / 'reference.tsv').symlink_to(reference)
        for destination, source in ((Path('/src/data'), staging), (Path('/src/outdir'), results)):
            require(not destination.exists() and not destination.is_symlink(),
                    'Image layout error: ' + str(destination) + ' already exists.')
            destination.symlink_to(source, target_is_directory=True)
        executable = Path('/src/CIBERSORTxFractions')
        require(executable.is_file() and os.access(str(executable), os.X_OK), 'Image error: CIBERSORTxFractions is not executable.')
        log('dimensions=genes:%d,cells:%d,cell_types:%d Start signature derivation.' % (len(input_genes), len(columns), len(counts)))
        command = [str(executable), '--username', args.username, '--token', args.token,
                   '--refsample', 'reference.tsv', '--single_cell', 'TRUE',
                   '--replicates', str(args.replicates), '--sampling', str(args.sampling),
                   '--fraction', str(args.fraction), '--G.min', str(args.min_genes),
                   '--G.max', str(args.max_genes), '--q.value', str(args.q_value),
                   '--k.max', str(args.max_condition_number), '--QN', 'FALSE', '--outdir', '/src/outdir']
        environment = dict(os.environ, OMP_NUM_THREADS=str(args.cpu), OPENBLAS_NUM_THREADS=str(args.cpu), MKL_NUM_THREADS=str(args.cpu))
        stopped = threading.Event()
        def heartbeat():
            while not stopped.wait(60):
                log('Signature derivation is active.')
        threading.Thread(target=heartbeat, daemon=True).start()
        try:
            process = subprocess.Popen(command, cwd='/src', env=environment, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1)
            for line in process.stdout:
                if not line.lstrip().startswith(('>[Options] token:', '>[Options] username:')):
                    emit(line.rstrip('\r\n'))
            status = process.wait()
        finally:
            stopped.set()
        log('Native command exited with status %d.' % status)
        require(status == 0, 'CIBERSORTxFractions failed with exit status %d.' % status)

        signature = one_output('*.bm.K*.txt')
        source = one_output('CIBERSORTx_cell_type_sourceGEP.txt')
        inferred = one_output('*_inferred_refsample.txt')
        classes = one_output('*_inferred_phenoclasses.txt')
        cell_types, signature_genes = matrix(signature, 'Output error: signature')
        source_types, source_genes = matrix(source, 'Output error: source GEPs')
        profiles, reference_genes = matrix(inferred, 'Output error: reference sample')
        require(set(cell_types) == set(counts), 'Output error: signature cell types do not match the reference.')
        require(source_types == cell_types, 'Output error: source GEP cell types do not match the signature.')
        require(signature_genes.issubset(source_genes), 'Output error: signature genes are absent from source GEPs.')
        require(len(profiles) == len(counts) * args.replicates, 'Output error: unexpected number of reference profiles.')
        with classes.open(newline='') as stream:
            class_rows = list(csv.reader(stream, delimiter='\t'))
        require(len(class_rows) == len(cell_types), 'Output error: unexpected number of phenotype classes.')
        require([row[0] for row in class_rows if row] == cell_types, 'Output error: phenotype class labels do not match the signature.')
        require(all(len(row) == len(profiles) + 1 and all(value in ('1', '2') for value in row[1:])
                    and row[1:].count('1') == args.replicates for row in class_rows),
                'Output error: invalid phenotype class codes or replicate counts.')
        require(all(sum(row[index] == '1' for row in class_rows) == 1 for index in range(1, len(profiles) + 1)),
                'Output error: each reference profile must have one phenotype class.')
        output_paths = {'signature_matrix': signature, 'source_geps': source,
                        'reference_sample': inferred, 'phenotype_classes': classes}
        for name, path in output_paths.items():
            shutil.copyfile(str(path), str(results / (name + '.tsv')))
        report = {'gene_count': len(signature_genes), 'cell_types': cell_types,
                  'source_gene_count': len(source_genes), 'reference_gene_count': len(reference_genes),
                  'reference_profile_count': len(profiles), 'native_files': {name: path.name for name, path in output_paths.items()},
                  'sha256': {name: digest(path) for name, path in output_paths.items()}}
        (root / 'output_validation.json').write_text(json.dumps(report, indent=2) + '\n')
        log('completion_time=%s dimensions=marker_genes:%d,cell_types:%d outputs=signature_matrix,source_geps,reference_sample,phenotype_classes,validation_reports status=0' % (now(), len(signature_genes), len(cell_types)))
    except (ValueError, OSError) as error:
        log('error=%s completion_time=%s status=failed' % (error, now()))
        sys.exit(1)
    finally:
        logfile.close()
    PY
  >>>

  output {
    File signature_matrix = "results/signature_matrix.tsv"
    File source_geps = "results/source_geps.tsv"
    File reference_sample = "results/reference_sample.tsv"
    File phenotype_classes = "results/phenotype_classes.tsv"
    Array[File] signature_heatmaps = glob("results/*.bm.K*.pdf")
    File input_validation = "input_validation.json"
    File output_validation = "output_validation.json"
    File run_log = "cibersortx_markers.log"
    File task_stdout = stdout()
    File task_stderr = stderr()
  }

  runtime {
    docker: "cibersortx/fractions@sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 0
    maxRetries: 0
  }
}
