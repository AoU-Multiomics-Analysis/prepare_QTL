"""Validated local files and explicit genomic coordinates."""
import csv
import gzip
import math
from pathlib import Path
import re

ASSOCIATION_FIELDS = ['variant_id', 'chrom', 'pos', 'phenotype_id', 'p_value', 'modality', 'dataset_id', 'cell_type']

def log(message):
    print('[regions] ' + message, flush=True)

def readable(value):
    value = str(value)
    if '://' in value:
        raise ValueError('localization error: unresolved URI: ' + value)
    path = Path(value).resolve()
    if not path.is_file():
        raise ValueError('input is not a readable local file: ' + value)
    with path.open('rb') as handle:
        handle.read(1)
    return path

def text_open(path, mode='rt'):
    return gzip.open(path, mode, newline='') if str(path).endswith('.gz') else open(path, mode, newline='')

def rows(path):
    path = readable(path)
    if path.suffix == '.parquet':
        import pyarrow.parquet as pq
        for batch in pq.ParquetFile(path).iter_batches(batch_size=65536):
            yield from batch.to_pylist()
    else:
        with text_open(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError('missing or duplicate table header: ' + str(path))
            for row in reader:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError('ragged table: ' + str(path))
                yield row

def write_table(path, fields, records):
    with text_open(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)

def lines(path):
    return readable(path).read_text().splitlines()

def local_file_list(path):
    values = lines(path)
    if any(not x for x in values):
        raise ValueError('empty path in file list')
    return [readable(x) for x in values]

def token(value, name):
    value = str(value)
    if not value or value in {'.', 'None', 'nan'} or re.search(r'\s', value):
        raise ValueError('invalid ' + name + ': ' + value)
    return value

def chromosome(value):
    value = str(value).removeprefix('chr')
    if value not in [str(i) for i in range(1,23)] + ['X']:
        raise ValueError('supported chromosomes are 1-22 and X; received ' + value)
    return value

def chrom_key(value):
    return 23 if value == 'X' else int(value)

def sizes(path):
    result = {}
    for line in lines(path):
        if not line or line.startswith('#'):
            continue
        cells = line.split()
        if cells[0].removeprefix('chr') not in [str(i) for i in range(1,23)] + ['X']:
            continue  # FASTA indices may also contain unplaced contigs, Y and MT.
        chrom = chromosome(cells[0])
        size = int(cells[1])
        if chrom in result or size <= 0:
            raise ValueError('duplicate chromosome or nonpositive chromosome size')
        result[chrom] = size
    if not result:
        raise ValueError('no supported chromosome sizes')
    return result

def assignments(path):
    result = {}
    for row in rows(path):
        iid = token(row['IID'], 'IID')
        group = token(row['ancestry'], 'ancestry')
        if iid in result:
            raise ValueError('duplicate ancestry assignment: ' + iid)
        result[iid] = group
    if not result:
        raise ValueError('empty ancestry assignments')
    return result

def seeds(path, chrom=None):
    result = {}
    for row in rows(path):
        if chrom is None or row['chrom'] == chrom:
            entry = (row['chrom'], int(row['pos']))
            if row['variant_id'] in result and result[row['variant_id']] != entry:
                raise ValueError('variant ID maps to multiple positions: ' + row['variant_id'])
            result[row['variant_id']] = entry
    return result

def interval(start, end, length, minimum):
    """Input and output are 0-based, half-open. Preserve all original positions."""
    start, end = max(0, start), min(length, end)
    target = min(minimum, length)
    if end - start < target:
        extra = target - (end - start)
        start -= extra // 2
        end += extra - extra // 2
        if start < 0:
            end -= start
            start = 0
        if end > length:
            start -= end - length
            end = length
    return start, end

def finite(value, name):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(name + ' must be finite')
    return number
