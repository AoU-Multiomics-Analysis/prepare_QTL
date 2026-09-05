#!/usr/bin/env python3
"""Pull-only smoke test of both methylation binaries on a small known cohort."""
import argparse
import csv
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='rust-image-smoke-') as temporary:
        directory = Path(temporary)
        samples = ['s1', 's2', 's3']
        for index, sample in enumerate(samples):
            (directory / (sample + '.bed')).write_text(
                '#chrom\tbegin\tend\tmod_score\ttype\tcov\n'
                f'chr1\t100\t101\t{20 + 30 * index}\tcombined\t20\n'
                f'chr1\t200\t201\t{30 + 20 * index}\tcombined\t30\n'
                'chrX\t300\t301\t50\tcombined\t20\n')
        (directory / 'manifest.tsv').write_text('sample_id\tfile_path\n' + ''.join(
            f'{sample}\t/data/{sample}.bed\n' for sample in samples))
        (directory / 'samples.tsv').write_text('sample_id\n' + '\n'.join(samples) + '\n')

        def run(binary, *arguments):
            subprocess.run(['docker', 'run', '--rm', '--user', f'{os.getuid()}:{os.getgid()}',
                            '--volume', f'{directory}:/data', args.image, binary, *arguments], check=True)

        run('methylation-filter', '--input-manifest', '/data/manifest.tsv',
            '--output-prefix', '/data/filter', '--num-threads', '1')
        with (directory / 'filter.methylation.sample_qc.tsv').open() as stream:
            qc = list(csv.DictReader(stream, delimiter='\t'))
        assert [row['sample_id'] for row in qc] == samples
        assert all(int(row['n_removed_by_chrom_filter']) == 1 for row in qc)
        calls = 'filter.methylation.autosome01.per_sample_qc.long.tsv.gz'
        # The cohort merger consumes one coordinate-sorted stream per sample.
        with gzip.open(directory / calls, 'rt') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            fields = reader.fieldnames
            call_rows = list(reader)
        assert len(call_rows) == 6
        for sample in samples:
            with gzip.open(directory / (sample + '.calls.tsv.gz'), 'wt') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
                writer.writeheader()
                writer.writerows(row for row in call_rows if row['sample_id'] == sample)
        (directory / 'calls.txt').write_text(''.join('/data/' + sample + '.calls.tsv.gz\n' for sample in samples))
        run('methylation-chromosome-merge', '--all-call-list', '/data/calls.txt',
            '--sample-qc', '/data/filter.methylation.sample_qc.tsv',
            '--cohort-samples', '/data/samples.tsv', '--total-samples', '3',
            '--chromosome', 'chr1', '--output-prefix', '/data/merge',
            '--min-methylation-mad', '0', '--skip-coverage-methylation-correlation')
        with gzip.open(directory / 'merge.methylation.raw.bed.gz', 'rt') as stream:
            rows = list(csv.reader(stream, delimiter='\t'))
        assert rows[0][4:] == samples
        assert len(rows) == 3, rows
        for row, start, expected in zip(rows[1:], ['100', '200'], [[.2, .5, .8], [.3, .5, .7]]):
            assert row[:3] == ['chr1', start, str(int(start) + 1)]
            assert all(math.isclose(float(value), wanted, abs_tol=1e-12)
                       for value, wanted in zip(row[4:], expected))
        with gzip.open(directory / 'merge.methylation.INT.bed.gz', 'rt') as stream:
            transformed = list(csv.reader(stream, delimiter='\t'))
        assert transformed[0] == rows[0] and len(transformed) == len(rows)
        assert all(all(math.isfinite(float(value)) for value in row[4:]) for row in transformed[1:])
    print('Both bundled Rust binaries passed the synthetic cohort check')


if __name__ == '__main__':
    main()
