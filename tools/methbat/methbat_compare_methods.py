#!/usr/bin/env python3
"""Compare MethBat output from the same donors before and after aggregation."""

import argparse
from collections import Counter
import csv
from decimal import Decimal
from pathlib import Path

OUTPUTS = ("meth_regions.bed", "combined_methyl.bedgraph", "asm.bedgraph")


def read_rows(path):
    with path.open(encoding="utf-8") as handle:
        return Counter(tuple(line.rstrip("\n").split("\t"))
                       for line in handle if line.strip() and not line.startswith("#"))


def compare(prefix, subset, chromosome, sample_count):
    with Path(f"{prefix}.method_comparison.tsv").open("w", encoding="utf-8", newline="") as report, \
            Path(f"{prefix}.method_differences.tsv").open("w", encoding="utf-8", newline="") as differences:
        summary = csv.writer(report, delimiter="\t", lineterminator="\n")
        detail = csv.writer(differences, delimiter="\t", lineterminator="\n")
        summary.writerow([
            "subset", "chromosome", "sample_count", "output", "joint_rows", "aggregate_rows",
            "exact_rows_match", "exact_matched_rows", "joint_only_rows", "aggregate_only_rows",
            "shared_coordinates", "max_abs_score_delta_pp",
        ])
        detail.writerow(["output", "method", "chromosome", "start", "end", "label_or_score"])
        for suffix in OUTPUTS:
            joint = read_rows(Path(f"{prefix}.{suffix}"))
            aggregate = read_rows(Path(f"{prefix}.aggregate.{suffix}"))
            joint_only, aggregate_only = joint - aggregate, aggregate - joint
            joint_coords = {row[:3]: row[3] for row in joint}
            aggregate_coords = {row[:3]: row[3] for row in aggregate}
            shared = joint_coords.keys() & aggregate_coords.keys()
            score_delta = "NA"
            if suffix.endswith(".bedgraph") and shared:
                score_delta = max(abs(Decimal(joint_coords[key]) - Decimal(aggregate_coords[key]))
                                  for key in shared)
            summary.writerow([
                subset, chromosome, sample_count, suffix, joint.total(), aggregate.total(),
                str(joint == aggregate).lower(), joint.total() - joint_only.total(),
                joint_only.total(), aggregate_only.total(),
                len(shared), score_delta,
            ])
            for method, rows in (("joint", joint_only), ("aggregate", aggregate_only)):
                for row in sorted(rows, key=lambda row: (row[0], int(row[1]), int(row[2]), row[3:])):
                    for _ in range(rows[row]):
                        detail.writerow([suffix, method, *row])
    print(f"Compared all three MethBat outputs for {subset} on {chromosome} with {sample_count} donors.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--subset", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--sample-count", type=int, required=True)
    args = parser.parse_args()
    compare(args.prefix, args.subset, args.chromosome, args.sample_count)


if __name__ == "__main__":
    main()
