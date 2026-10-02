version 1.0

workflow CIBERSORTxFractions {
  input {
    File mixture
    File signature
    String username
    String token
    File? refsample
    File? source_geps
    Boolean smode = false
    Boolean quantile_normalization = false
    Int permutations = 100
    Int threads = 8
    Int memory_gb = 16
    Int disk_gb = 30
  }

  call RunFractions {
    input:
      mixture = mixture,
      signature = signature,
      username = username,
      token = token,
      refsample = refsample,
      source_geps = source_geps,
      smode = smode,
      quantile_normalization = quantile_normalization,
      permutations = permutations,
      threads = threads,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File fractions = RunFractions.fractions
    File fractions_only = RunFractions.fractions_only
    File mixture_for_hires = RunFractions.mixture_for_hires
    File signature_for_hires = RunFractions.signature_for_hires
    Array[File] adjusted_mixture = RunFractions.adjusted_mixture
    Array[File] adjusted_signature = RunFractions.adjusted_signature
    File input_validation = RunFractions.input_validation
    File output_validation = RunFractions.output_validation
    File run_log = RunFractions.run_log
    File task_stdout = RunFractions.task_stdout
    File task_stderr = RunFractions.task_stderr
  }
}

task RunFractions {
  input {
    File mixture
    File signature
    String username
    String token
    File? refsample
    File? source_geps
    Boolean smode = false
    Boolean quantile_normalization = false
    Int permutations = 100
    Int threads = 8
    Int memory_gb = 16
    Int disk_gb = 30
  }

  command <<<
    set -euo pipefail
    # File values stay typed until the command is rendered. Optional flags are
    # included only when the corresponding File input is present.
    python3 - \
      --mixture '~{sub(mixture, "'", "'\"'\"'")}' \
      --signature '~{sub(signature, "'", "'\"'\"'")}' \
      --username '~{sub(username, "'", "'\"'\"'")}' \
      --token '~{sub(token, "'", "'\"'\"'")}' \
      ~{if defined(refsample) then "--refsample '" + sub(select_first([refsample]), "'", "'\"'\"'") + "'" else ""} \
      ~{if defined(source_geps) then "--source-geps '" + sub(select_first([source_geps]), "'", "'\"'\"'") + "'" else ""} \
      --smode ~{if smode then "TRUE" else "FALSE"} \
      --qn ~{if quantile_normalization then "TRUE" else "FALSE"} \
      --permutations ~{permutations} --threads ~{threads} \
      --memory-gb ~{memory_gb} --disk-gb ~{disk_gb} <<'PY'
    import argparse, csv, datetime, gzip, hashlib, json, math, os, shutil, subprocess, sys, tempfile, threading, zlib
    from pathlib import Path

    parser = argparse.ArgumentParser()
    for name in ('mixture', 'signature', 'username', 'token'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--refsample')
    parser.add_argument('--source-geps')
    for name in ('smode', 'qn'):
        parser.add_argument('--' + name, choices=['TRUE', 'FALSE'], required=True)
    for name in ('permutations', 'threads', 'memory-gb', 'disk-gb'):
        parser.add_argument('--' + name, type=int, required=True)
    args = parser.parse_args()
    root = Path.cwd()
    logfile = (root / 'fractions.log').open('w', buffering=1)
    def emit(message):
        print(message, flush=True)
        logfile.write(message + '\n')
    def timestamp():
        return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    def log(message):
        emit('[%s] stage=CIBERSORTxFractions %s' % (timestamp(), message))
    def require(condition, message):
        if not condition:
            raise SystemExit(message)
    def readable(value, label):
        require(value is not None and '://' not in value,
                'Localization error: ' + label + ' must be a readable localized File.')
        path = Path(value).resolve()
        require(path.is_file() and os.access(str(path), os.R_OK),
                'Localization error: ' + label + ' must be a readable localized File.')
        return path
    def table(path, context, allow_negative=False, unique_columns=True):
        with path.open(newline='') as stream:
            reader = csv.reader(stream, delimiter='\t')
            header = next(reader, [])
            require(len(header) > 1, context + ': matrix has no numeric columns.')
            columns = header[1:]
            require(all(x and x == x.strip() for x in columns), context + ': empty or padded column label.')
            if unique_columns:
                require(len(set(columns)) == len(columns), context + ': duplicate column labels.')
            rows = []
            seen = set()
            negative = 0
            totals = [0.0] * len(columns)
            for row in reader:
                require(len(row) == len(header), context + ': matrix is not rectangular.')
                require(row[0] and row[0] == row[0].strip() and row[0] not in seen,
                        context + ': empty, padded, or duplicate gene label.')
                seen.add(row[0])
                try:
                    values = [float(x) for x in row[1:]]
                except ValueError:
                    raise SystemExit(context + ': matrix contains nonnumeric expression.')
                require(all(math.isfinite(x) for x in values), context + ': nonfinite expression.')
                negative += sum(x < 0 for x in values)
                require(allow_negative or all(x >= 0 for x in values), context + ': negative expression; supply unadjusted linear expression.')
                totals = [a + b for a, b in zip(totals, values)]
                rows.append(row)
            require(bool(rows), context + ': matrix has no gene rows.')
            return header, rows, seen, negative, totals
    def write_table(path, header, rows):
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(header)
            writer.writerows(rows)
    def digest(path):
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(chunk)
        return h.hexdigest()

    log('start_time=%s Check localized inputs.' % timestamp())
    require(args.permutations >= 0, 'Input error: permutations must be nonnegative.')
    require(all(getattr(args, x) > 0 for x in ('threads', 'memory_gb', 'disk_gb')),
            'Input error: CPU, memory, and disk requests must be positive.')
    require(bool(args.username.strip()) and bool(args.token.strip()), 'Input error: username and token are required.')
    smode = args.smode == 'TRUE'
    require(not smode or (args.refsample and args.source_geps),
            'Input error: S-mode requires refsample and source_geps with a prebuilt signature.')
    paths = {name: readable(getattr(args, name), name) for name in ('mixture', 'signature')}
    for name in ('refsample', 'source_geps'):
        value = getattr(args, name)
        if value is not None:
            paths[name] = readable(value, name)
    supplied_paths = dict(paths)
    with paths['mixture'].open('rb') as stream:
        mixture_is_gzip = stream.read(2) == b'\x1f\x8b'
    if mixture_is_gzip:
        log('Decompress the localized gzip mixture to task disk.')
        uncompressed = Path(tempfile.mkdtemp(prefix='uncompressed_mixture_', dir=str(root))) / 'mixture.txt'
        try:
            with gzip.open(str(paths['mixture']), 'rb') as source, uncompressed.open('wb') as destination:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
        except (OSError, EOFError, zlib.error) as error:
            message = 'Input error: mixture gzip could not be decompressed: ' + str(error)
            log(message)
            raise SystemExit(message)
        paths['mixture'] = uncompressed
        log('Mixture decompression completed; check the uncompressed matrix.')
    mh, mr, mixture_genes, _, mt = table(paths['mixture'], 'Input error')
    sh, sr, signature_genes, _, st = table(paths['signature'], 'Input error')
    samples, cells = mh[1:], sh[1:]
    require(all(x > 0 for x in mt), 'Input error: each mixture sample must have positive total expression.')
    require(all(x > 0 for x in st), 'Input error: each signature cell type must have positive total expression.')
    require(bool(mixture_genes & signature_genes), 'Input error: mixture and signature share no genes.')
    for name in ('refsample', 'source_geps'):
        if name in paths:
            h, rows, genes, _, totals = table(paths[name], 'Input error', unique_columns=(name != 'refsample'))
            require(bool(genes & mixture_genes), 'Input error: ' + name + ' and mixture share no genes.')
            if name == 'source_geps':
                require(h[1:] == cells, 'Input error: source_geps cell labels and order must match the signature.')

    results, inputs = root / 'results', root / 'inputs'
    results.mkdir()
    inputs.mkdir()
    for name, path in paths.items():
        (inputs / (name + '.txt')).symlink_to(path)
    for destination, source in ((Path('/src/data'), inputs), (Path('/src/outdir'), results)):
        if destination.is_dir() and not destination.is_symlink() and not any(destination.iterdir()):
            destination.rmdir()
        require(not destination.exists() and not destination.is_symlink(), 'Staging error: ' + str(destination) + ' already contains files.')
        destination.symlink_to(source, target_is_directory=True)
    input_report = {'samples': samples, 'sample_count': len(samples), 'cell_types': cells,
                    'mixture_gene_count': len(mixture_genes), 'signature_gene_count': len(signature_genes),
                    'shared_gene_count': len(mixture_genes & signature_genes), 'smode': smode,
                    'quantile_normalization': args.qn == 'TRUE', 'permutations': args.permutations,
                    'sha256': {name: digest(path) for name, path in supplied_paths.items()},
                    'mixture_compression': 'gzip' if mixture_is_gzip else 'none',
                    'mixture_bytes': supplied_paths['mixture'].stat().st_size,
                    'uncompressed_mixture_bytes': paths['mixture'].stat().st_size,
                    'uncompressed_mixture_sha256': digest(paths['mixture'])}
    (root / 'input_validation.json').write_text(json.dumps(input_report, indent=2) + '\n')
    log('dimensions=mixture_genes:%d,samples:%d,signature_genes:%d,cell_types:%d' %
        (len(mixture_genes), len(samples), len(signature_genes), len(cells)))

    command = ['./CIBERSORTxFractions', '--username', args.username, '--token', args.token,
               '--mixture', 'mixture.txt', '--sigmatrix', 'signature.txt',
               '--QN', args.qn, '--perm', str(args.permutations), '--absolute', 'FALSE',
               '--rmbatchSmode', args.smode, '--outdir', '/src/outdir']
    for name, flag in (('refsample', '--refsample'), ('source_geps', '--sourceGEPs')):
        if name in paths:
            command.extend([flag, name + '.txt'])
    environment = dict(os.environ)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        environment[key] = str(args.threads)
    require(os.access('/src/CIBERSORTxFractions', os.X_OK), 'Image error: CIBERSORTxFractions is not executable.')
    log('Start relative fraction estimation: %d samples; %d cell types; permutations=%d; S-mode=%s.' %
        (len(samples), len(cells), args.permutations, args.smode))
    # The native Fractions interface has no --threads option. Limit threaded
    # libraries through the environment; this does not guarantee parallel fits.
    secrets = [part for value in (args.username, args.token) for part in value.splitlines() if part]
    finished = threading.Event()
    def heartbeat():
        minutes = 0
        while not finished.wait(60):
            minutes += 1
            log('Native command active for %d minutes.' % minutes)
    monitor = threading.Thread(target=heartbeat, daemon=True)
    monitor.start()
    try:
        process = subprocess.Popen(command, cwd='/src', env=environment, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1)
        for line in process.stdout:
            if line.startswith(('>[Options] token:', '>[Options] username:')):
                continue
            if any(secret in line for secret in secrets):
                emit('[credential text omitted]')
            else:
                emit(line.rstrip('\n'))
        status = process.wait()
    finally:
        finished.set()
        monitor.join()
    log('Native command exited with status %d.' % status)
    if status != 0:
        log('completion_time=%s native_exit_status=%d' % (timestamp(), status))
        sys.exit(status if status > 0 else 128 - status)

    log('Check fraction outputs and prepare HiRes inputs.')
    fraction_path = results / ('CIBERSORTx_Adjusted.txt' if smode else 'CIBERSORTx_Results.txt')
    require(fraction_path.is_file(), 'Output error: expected fraction table is missing.')
    with fraction_path.open(newline='') as stream:
        rows = list(csv.reader(stream, delimiter='\t'))
    require(len(rows) > 1 and len(set(rows[0])) == len(rows[0]), 'Output error: empty or duplicate fraction columns.')
    fh = rows[0]
    require(set(fh[1:]) == set(cells) | {'P-value', 'Correlation', 'RMSE'} or set(fh[1:]) == set(cells),
            'Output error: fraction cell labels do not match the signature.')
    ids = [row[0] for row in rows[1:] if row]
    require(len(ids) == len(samples) and len(set(ids)) == len(ids) and set(ids) == set(samples),
            'Output error: fraction sample IDs do not match the mixture.')
    by_sample = {}
    indexes = [fh.index(cell) for cell in cells]
    for row in rows[1:]:
        require(len(row) == len(fh), 'Output error: fraction table is not rectangular.')
        try:
            fractions = [float(row[i]) for i in indexes]
        except ValueError:
            raise SystemExit('Output error: nonnumeric fraction.')
        require(all(math.isfinite(x) and 0 <= x <= 1 for x in fractions), 'Output error: nonfinite or out-of-range fraction.')
        require(abs(sum(fractions) - 1) <= 1e-4, 'Output error: relative fractions do not sum to one.')
        by_sample[row[0]] = [row[0]] + [row[i] for i in indexes]
    write_table(results / 'fractions_only.txt', ['Mixture'] + cells, [by_sample[x] for x in samples])
    shutil.copyfile(fraction_path, results / 'fractions.txt')

    export_mixture, export_signature = paths['mixture'], paths['signature']
    if smode:
        export_mixture = results / 'CIBERSORTx_Mixtures_Adjusted.txt'
        # Native versions can name this after the staged signature filename.
        candidates = [x for x in results.glob('*_Adjusted.txt')
                      if x.name not in ('CIBERSORTx_Adjusted.txt', 'CIBERSORTx_Mixtures_Adjusted.txt')]
        require(export_mixture.is_file() and len(candidates) == 1,
                'Output error: S-mode adjusted mixture or signature is missing or ambiguous.')
        export_signature = candidates[0]
    eh, er, eg, en, et = table(export_mixture, 'Output error', allow_negative=smode)
    ah, ar, ag, an, at = table(export_signature, 'Output error', allow_negative=smode)
    require(eh[1:] == samples and ah[1:] == cells, 'Output error: adjusted sample or cell labels changed.')
    shared_rows = [row for row in ar if row[0] in eg]
    require(bool(shared_rows), 'Output error: exported signature and mixture share no genes.')
    write_table(results / 'signature_shared_genes.txt', ah, shared_rows)
    if smode:
        shutil.copyfile(export_mixture, results / 'mixture_adjusted.txt')
        shutil.copyfile(export_signature, results / 'signature_adjusted.txt')
    # Always collect plain text for the downstream HiRes task, including when
    # the unadjusted input was compressed. The source input stays unchanged.
    shutil.copyfile(export_mixture, results / 'mixture_for_hires.txt')
    output_report = {'sample_count': len(samples), 'samples': samples, 'cell_types': cells,
                     'fraction_table': fraction_path.name, 'fraction_row_sum_tolerance': 1e-4,
                     'signature_shared_genes': len(shared_rows),
                     'signature_removed_genes': len(ar) - len(shared_rows),
                     'negative_adjusted_mixture_values': en, 'negative_adjusted_signature_values': an,
                     'smode': smode, 'fraction_scale': 'relative; rows sum to one'}
    (root / 'output_validation.json').write_text(json.dumps(output_report, indent=2) + '\n')
    log('outputs=%s' % ','.join(sorted(path.name for path in results.iterdir() if path.is_file())))
    log('Fraction task complete: %d samples; %d cell types; %d shared signature genes.' %
        (len(samples), len(cells), len(shared_rows)))
    log('completion_time=%s exit_status=0' % timestamp())
    PY
  >>>

  output {
    File fractions = "results/fractions.txt"
    File fractions_only = "results/fractions_only.txt"
    File mixture_for_hires = "results/mixture_for_hires.txt"
    File signature_for_hires = "results/signature_shared_genes.txt"
    Array[File] adjusted_mixture = glob("results/mixture_adjusted.txt")
    Array[File] adjusted_signature = glob("results/signature_adjusted.txt")
    File input_validation = "input_validation.json"
    File output_validation = "output_validation.json"
    File run_log = "fractions.log"
    File task_stdout = stdout()
    File task_stderr = stderr()
  }

  runtime {
    docker: "cibersortx/fractions@sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc"
    cpu: threads
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} SSD"
    bootDiskSizeGb: 20
    preemptible: 0
    maxRetries: 0
  }
}
