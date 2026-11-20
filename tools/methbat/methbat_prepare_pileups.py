#!/usr/bin/env python3
"""Convert staged legacy CpG tracks and make a one-row MethBat cohort."""

import argparse
from contextlib import contextmanager
import csv
import gzip
from math import modf
from heapq import merge
from itertools import groupby, islice
from pathlib import Path
import re
import subprocess

TRACKS = ("Total", "hap1", "hap2")
LEGACY_TRACKS = ("combined", "hap1", "hap2")
MANIFEST_FIELDS = ("identifier", "filename", "labels")
BED_HEADER = (
    "##base_modification=5mC\n"
    "#chrom\tstart\tend\tname\tscore\tstrand\tmod_score\ttype\tcov\tmod_count\tunmod_count"
    "\tinferred_unmod_count\tdiff_base_count\tavg_mod_score\tavg_unmod_score\n"
)
# MethBat v1.2.0 loads fractions as rounded integers on this scale.
FRACTION_SCALE = 10_000


def round_ratio(numerator, denominator):
    quotient, remainder = divmod(numerator, denominator)
    return quotient + int(2 * remainder >= denominator)


@contextmanager
def bgzip_writer(path):
    with path.open("wb") as compressed:
        with subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE,
                              stdout=compressed, text=True, encoding="utf-8") as bgzip:
            with bgzip.stdin as handle:
                handle.write(BED_HEADER)
                yield handle
            if bgzip.wait() != 0:
                raise ValueError(f"bgzip could not write {path}.")


def legacy_rows(path, chromosome, track):
    previous = -1
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            try:
                if len(columns) not in (9, 10) or columns[0] != chromosome:
                    raise ValueError
                start, end = map(int, columns[1:3])
                coverage, modified, unmodified = map(int, columns[5:8])
                if (start <= previous or end != start + 1 or not 0 <= coverage < 2**32
                        or min(modified, unmodified) < 0 or modified + unmodified != coverage):
                    raise ValueError
            except ValueError as error:
                raise ValueError(f"Invalid or unsorted pb-CpG-tools row in {path}, line {number}.") from error
            previous = start
            if coverage:
                # Model-mode files supply estimated counts and no mean probabilities.
                averages = columns[8:10] if len(columns) == 10 else ("0.000", "0.000")
                yield start, track, coverage, modified, averages


def bed_row(chromosome, start, track, coverage, modified, averages=("0.000", "0.000")):
    score = round_ratio(1000 * modified, coverage)
    percent = f"{score // 10}.{score % 10}"
    return (f"{chromosome}\t{start}\t{start + 1}\tm\t{score}\t.\t{percent}\t{TRACKS[track]}"
            f"\t{coverage}\t{modified}\t{coverage - modified}\t0\t0\t{averages[0]}\t{averages[1]}\n")


def convert_sample(prefix, chromosome):
    output = Path(f"{prefix}.5mC.bed.gz")
    sources = [legacy_rows(Path(f"{prefix}.{name}.bed.gz"), chromosome, track)
               for track, name in enumerate(LEGACY_TRACKS)]
    with bgzip_writer(output) as handle:
        for start, track, coverage, modified, averages in merge(*sources):
            handle.write(bed_row(chromosome, start, track, coverage, modified, averages))
    return output


def converted_sites(path, *, exact_ratio=False):
    def rows(handle):
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.split("\t", 11)
            modified, coverage = int(columns[9]), int(columns[8])
            if exact_ratio:
                # Match joint-segment's floating-point conversion of each donor.
                remainder, whole = modf((modified / coverage) * FRACTION_SCALE)
                fraction = int(whole) + int(remainder >= 0.5)
            else:
                fraction = round_ratio(FRACTION_SCALE * modified, coverage)
            yield int(columns[1]), TRACKS.index(columns[7]), fraction

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for start, calls in groupby(rows(handle), key=lambda row: row[0]):
            fractions = [None, None, None]
            for _, track, value in calls:
                fractions[track] = value
            yield start, fractions


def write_aggregate(pileups, chromosome, output, *, exact_ratio=False):
    def aggregate_row(start, track, total, count):
        if exact_ratio:
            # The count columns hold the sum of donor fractions and its scale.
            return bed_row(chromosome, start, track, FRACTION_SCALE * count, total)
        return bed_row(chromosome, start, track, FRACTION_SCALE, round_ratio(total, count))

    candidate_sites = combined_sites = asm_sites = 0
    n_donors = len(pileups)
    merged = merge(*(converted_sites(path, exact_ratio=exact_ratio) for path in pileups),
                   key=lambda site: site[0])
    with bgzip_writer(output) as handle:
        for start, observations in groupby(merged, key=lambda site: site[0]):
            candidate_sites += 1
            combined_count = combined_sum = paired_count = low_sum = high_sum = 0
            for _, (total, hap1, hap2) in observations:
                if total is not None:
                    combined_count += 1
                    combined_sum += total
                if hap1 is not None and hap2 is not None:
                    paired_count += 1
                    low, high = sorted((hap1, hap2))
                    low_sum += low
                    high_sum += high
            # Exact count thresholds for --min-sample-frac 0.9 and --min-asm-frac 0.6.
            if 10 * combined_count < 9 * n_donors:
                continue
            combined_sites += 1
            handle.write(aggregate_row(start, 0, combined_sum, combined_count))
            if 10 * paired_count >= 6 * n_donors:
                asm_sites += 1
                for track, total in ((1, low_sum), (2, high_sum)):
                    handle.write(aggregate_row(start, track, total, paired_count))
    return candidate_sites, combined_sites, asm_sites


def prepare(args):
    with args.collection.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        if tuple(next(reader, [])) != MANIFEST_FIELDS:
            raise ValueError("Invalid chromosome bundle collection header.")
        samples = list(islice(reader, args.sample_limit) if args.sample_limit else reader)
    if not samples or (args.sample_limit and len(samples) != args.sample_limit):
        raise ValueError("The chromosome bundle does not contain enough samples.")

    donor_collection = Path("selected.tsv")
    pileups = []
    with donor_collection.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(MANIFEST_FIELDS) + "\n")
        for identifier, prefix, label in samples:
            if not re.fullmatch(r"tracks/sample-[0-9]+", prefix):
                raise ValueError(f"Invalid staged track prefix {prefix} in the chromosome bundle.")
            pileup = convert_sample(prefix, args.chromosome)
            pileups.append(pileup)
            handle.write(f"{identifier}\t{pileup}\t{label}\n")

    aggregate = Path(f"{args.output_prefix}.aggregate.5mC.bed.gz")
    candidates, combined, asm = write_aggregate(pileups, args.chromosome, aggregate,
                                                 exact_ratio=args.exact_ratio)
    if not args.exact_ratio:
        Path("aggregate.tsv").write_text(
            "\t".join(MANIFEST_FIELDS) + f"\ncohort\t{aggregate}\t\n", encoding="utf-8"
        )
    Path("sample_count.txt").write_text(f"{len(samples)}\n", encoding="utf-8")
    with Path(f"{args.output_prefix}.aggregation_qc.tsv").open("w", encoding="utf-8") as handle:
        handle.write("chromosome\tsample_count\tcandidate_sites\tcombined_pass_sites\tpaired_haplotype_pass_sites\n")
        handle.write(f"{args.chromosome}\t{len(samples)}\t{candidates}\t{combined}\t{asm}\n")
    print(f"Prepared {len(samples)} donor pileups. Retained {combined} of {candidates} CpGs, with {asm} passing the paired-haplotype filter.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--sample-limit", type=int, required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--exact-ratio", action="store_true",
                        help="Preserve cohort ratios for methbat segment.")
    args = parser.parse_args()
    try:
        prepare(args)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Pileup preparation error. {error}\n")


if __name__ == "__main__":
    main()
