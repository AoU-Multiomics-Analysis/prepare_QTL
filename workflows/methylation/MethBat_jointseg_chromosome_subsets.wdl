version 1.0

workflow MethBat_jointseg_chromosome_subsets {
    meta {
        description: "Run MethBat v0.17 joint segmentation on disjoint cohorts and measure boundary shifts by sample size and chromosome."
    }

    input {
        File manifest_file
        Array[String]+ chromosomes = ["chr21"]
        Int subset_count = 2
        Array[Int] subset_sizes = if subset_count == 1 then [] else [25, 50, 100]
        Int split_seed = 42
        Float min_reciprocal_overlap = 0.5

        String? gcs_billing_project
        Int transfer_workers = 4
        Int staging_disk_gb = 100
        Int threads = 4
        Int mem_gb = 32
        Int disk_gb = 50
        Boolean condense_bed_labels = false
        Boolean enable_nodata_segments = false
        String staging_docker = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        String methbat_docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v2"
    }

    parameter_meta {
        manifest_file: "TSV with identifier, gs:// pb-CpG-tools prefix, and optional labels. Each prefix needs combined, hap1, and hap2 BEDs and their .tbi indexes. Use immutable source paths for Terra call caching."
        chromosomes: "Exact contig names from the input BED indexes."
        subset_count: "Number of disjoint label-balanced groups. Use 1 for the full cohort."
        subset_sizes: "Number of samples to take from each group. Sizes are nested within a group. Use an empty array to run each entire group."
        split_seed: "Random seed for reproducible subset membership."
        min_reciprocal_overlap: "Minimum fraction of each region covered by a same-label match in the other subset."
    }

    call SplitMethBatCohort {
        input:
            manifest_file = manifest_file,
            subset_count = subset_count,
            subset_sizes = subset_sizes,
            split_seed = split_seed,
            docker_image = staging_docker
    }

    Array[Int] sample_limits = if length(subset_sizes) > 0 then subset_sizes else [0]

    scatter (subset_manifest in SplitMethBatCohort.subset_manifests) {
        String subset_name = basename(subset_manifest, ".tsv")
        call StageMethBatChromosomes {
            input:
                subset_manifest = subset_manifest,
                chromosomes = chromosomes,
                transfer_workers = transfer_workers,
                disk_gb = staging_disk_gb,
                gcs_billing_project = gcs_billing_project,
                docker_image = staging_docker
        }

        scatter (run in cross(StageMethBatChromosomes.chromosome_bundles, sample_limits)) {
            call SegmentMethBatChromosome {
                input:
                    chromosome_bundle = run.left,
                    subset_name = subset_name,
                    chromosome = basename(run.left, ".tar"),
                    sample_limit = run.right,
                    threads = threads,
                    mem_gb = mem_gb,
                    disk_gb = disk_gb,
                    condense_bed_labels = condense_bed_labels,
                    enable_nodata_segments = enable_nodata_segments,
                    docker_image = methbat_docker
            }
        }
    }

    call CompareMethBatBoundaries {
        input:
            region_beds = flatten(SegmentMethBatChromosome.meth_regions_bed),
            run_keys = flatten(SegmentMethBatChromosome.run_key),
            min_reciprocal_overlap = min_reciprocal_overlap,
            docker_image = staging_docker
    }

    output {
        File membership_tsv = SplitMethBatCohort.membership_tsv
        Array[File] subset_manifests = SplitMethBatCohort.subset_manifests
        # Arrays are indexed first by group and then by chromosome/size run.
        Array[String] subsets = subset_name
        Array[Array[String]] result_chromosomes = SegmentMethBatChromosome.result_chromosome
        Array[Array[Int]] result_sample_counts = SegmentMethBatChromosome.sample_count
        Array[Array[File]] meth_regions_bed = SegmentMethBatChromosome.meth_regions_bed
        Array[Array[File]] combined_methyl_bedgraph = SegmentMethBatChromosome.combined_methyl_bedgraph
        Array[Array[File]] asm_bedgraph = SegmentMethBatChromosome.asm_bedgraph
        File boundary_summary_tsv = CompareMethBatBoundaries.boundary_summary
        File boundary_shifts_tsv = CompareMethBatBoundaries.boundary_shifts
    }
}

task SplitMethBatCohort {
    input {
        File manifest_file
        Int subset_count
        Array[Int] subset_sizes = []
        Int split_seed
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
    }

    command <<<
        set -euo pipefail

        manifest='~{sub("~{manifest_file}", "'", "'\"'\"'")}'
        if [[ "$manifest" == *://* || ! -r "$manifest" ]]; then
            echo "Manifest localization error. Cannot read the local input manifest." >&2
            exit 1
        fi

        python3 - --manifest "$manifest" --subsets ~{subset_count} --seed ~{split_seed} \
            --sizes ~{sep=" " subset_sizes} <<'PY'
        import argparse
        from collections import defaultdict
        from pathlib import Path
        import random
        import re

        parser = argparse.ArgumentParser()
        parser.add_argument("--manifest", type=Path, required=True)
        parser.add_argument("--subsets", type=int, required=True)
        parser.add_argument("--seed", type=int, required=True)
        parser.add_argument("--sizes", type=int, nargs="*", default=[])
        args = parser.parse_args()

        samples_by_label = defaultdict(list)
        identifiers = set()
        prefixes = set()
        first_row = True
        with args.manifest.open(encoding="utf-8-sig") as manifest:
            for line_number, line in enumerate(manifest, start=1):
                if not line.strip():
                    continue
                columns = line.rstrip("\n").split("\t")
                if first_row and columns in (["identifier", "filename"], ["identifier", "filename", "labels"]):
                    first_row = False
                    continue
                first_row = False
                if len(columns) not in (2, 3):
                    raise SystemExit(f"Manifest line {line_number} must have two or three TSV fields.")
                identifier, prefix = columns[:2]
                label = columns[2] if len(columns) == 3 else ""
                if not identifier.strip() or "\0" in identifier or "\0" in label:
                    raise SystemExit(f"Invalid identifier or label on manifest line {line_number}.")
                # Reject GCS glob syntax and paths that end in a filename suffix.
                if not re.fullmatch(r"gs://[^/*?\[\]#\s\x00]+/[^*?\[\]#\s\x00]+", prefix) or prefix.endswith(
                    ("/", ".bed", ".bed.gz", ".tbi", ".csi")
                ):
                    raise SystemExit(f"Manifest line {line_number} requires a gs:// pb-CpG-tools prefix.")
                if identifier in identifiers or prefix in prefixes:
                    raise SystemExit(f"Duplicate sample identifier or prefix on manifest line {line_number}.")
                identifiers.add(identifier)
                prefixes.add(prefix)
                samples_by_label[label].append((identifier, prefix, label))

        if not 1 <= args.subsets <= len(identifiers):
            raise SystemExit("Subset count must be between 1 and the number of samples.")
        if any(size < 1 for size in args.sizes) or len(args.sizes) != len(set(args.sizes)):
            raise SystemExit("Subset sizes must be unique positive integers.")
        max_size = max(args.sizes, default=0)
        needed_samples = args.subsets * max_size
        if needed_samples > len(identifiers):
            raise SystemExit(f"The requested disjoint groups require {needed_samples} unique samples. "
                             f"The manifest contains {len(identifiers)}.")

        groups = [defaultdict(list) for _ in range(args.subsets)]
        rng = random.Random(args.seed)
        offset = 0
        for label in sorted(samples_by_label):
            samples = sorted(samples_by_label[label])
            rng.shuffle(samples)
            for index, sample in enumerate(samples):
                groups[(index + offset) % args.subsets][label].append(sample)
            offset = (offset + len(samples)) % args.subsets

        output = Path("subsets")
        output.mkdir()
        with Path("membership.tsv").open("w", encoding="utf-8") as membership:
            membership.write("identifier\tlabels\tsubset\tselection_rank\n")
            for number, group in enumerate(groups, start=1):
                subset = f"subset-{number:06d}"
                if max_size:
                    # Interleave labels by their rank within the group so early prefixes
                    # retain approximately the same label mix as the full group.
                    ranked = sorted(
                        ((index + 0.5) / len(samples), label, sample)
                        for label, samples in group.items()
                        for index, sample in enumerate(samples)
                    )
                    selected = [sample for _, _, sample in ranked[:max_size]]
                else:
                    selected = sorted(sample for samples in group.values() for sample in samples)
                with (output / f"{subset}.tsv").open("w", encoding="utf-8") as handle:
                    handle.write("identifier\tfilename\tlabels\n")
                    for rank, (identifier, prefix, label) in enumerate(selected, start=1):
                        handle.write(f"{identifier}\t{prefix}\t{label}\n")
                        membership.write(f"{identifier}\t{label}\t{subset}\t{rank}\n")
                print(f"Wrote {subset} with {len(selected)} samples.", flush=True)
        PY
    >>>

    runtime {
        docker: docker_image
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
    }

    output {
        Array[File]+ subset_manifests = glob("subsets/subset-*.tsv")
        File membership_tsv = "membership.tsv"
    }
}

task StageMethBatChromosomes {
    input {
        File subset_manifest
        Array[String]+ chromosomes
        Int transfer_workers
        Int disk_gb = 100
        String? gcs_billing_project
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
    }

    command <<<
        set -euo pipefail

        manifest='~{sub("~{subset_manifest}", "'", "'\"'\"'")}'
        # Keep the generated path as a File until Cromwell renders the command.
        chromosome_list="~{write_lines(chromosomes)}"
        if [[ "$manifest" == *://* || ! -r "$manifest" ]]; then
            echo "Staging localization error. Cannot read the local subset manifest." >&2
            exit 1
        fi
        if [[ "$chromosome_list" == *://* || ! -r "$chromosome_list" ]]; then
            echo "Staging localization error. Cannot read the generated chromosome list." >&2
            exit 1
        fi

        for tool in gsutil tabix bgzip; do
            command -v "$tool" >/dev/null || { echo "Staging image lacks $tool." >&2; exit 1; }
        done
        billing_project='~{sub(select_first([gcs_billing_project, ""]), "'", "'\"'\"'")}'
        billing_args=()
        if [[ -n "$billing_project" ]]; then
            billing_args+=(--billing-project "$billing_project")
        fi

        python3 - --manifest "$manifest" --chromosome-list "$chromosome_list" \
            --chromosome-count ~{length(chromosomes)} --workers ~{transfer_workers} \
            "${billing_args[@]}" <<'PY'
        import argparse
        from concurrent.futures import ThreadPoolExecutor
        from pathlib import Path
        import re
        import subprocess
        import tarfile
        import tempfile

        parser = argparse.ArgumentParser()
        parser.add_argument("--manifest", type=Path, required=True)
        parser.add_argument("--chromosome-list", type=Path, required=True)
        parser.add_argument("--chromosome-count", type=int, required=True)
        parser.add_argument("--workers", type=int, required=True)
        parser.add_argument("--billing-project")
        args = parser.parse_args()

        chromosomes = args.chromosome_list.read_text(encoding="utf-8").splitlines()
        if len(chromosomes) != args.chromosome_count or any(
            not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", chrom) for chrom in chromosomes
        ):
            raise SystemExit("Invalid chromosome name. Use the exact contig name from the BED index.")
        if len(chromosomes) != len(set(chromosomes)):
            raise SystemExit("Chromosome names must be unique.")
        if args.workers < 1:
            raise SystemExit("transfer_workers must be at least 1.")

        with args.manifest.open(encoding="utf-8") as handle:
            if handle.readline().rstrip("\n") != "identifier\tfilename\tlabels":
                raise SystemExit("Invalid subset manifest header.")
            samples = [line.rstrip("\n").split("\t") for line in handle]
        if not samples:
            raise SystemExit("Subset manifest is empty.")

        tracks = ("combined", "hap1", "hap2")
        slice_root = Path("slices")
        bundles = Path("bundles")
        bundles.mkdir()
        for chrom in chromosomes:
            (slice_root / chrom / "tracks").mkdir(parents=True)

        gsutil = ["gsutil", "-q"]
        if args.billing_project:
            gsutil += ["-u", args.billing_project]

        def stage_sample(item):
            index, (identifier, prefix, label) = item
            local_stem = f"sample-{index:06d}"
            covered = set()
            with tempfile.TemporaryDirectory(prefix="download-", dir=".") as directory:
                temporary = Path(directory)
                for track in tracks:
                    # Download a full track once, even when several chromosomes are requested.
                    source = f"{prefix}.{track}.bed.gz"
                    local = temporary / source.rsplit("/", 1)[1]
                    subprocess.run([*gsutil, "cp", source, source + ".tbi", directory + "/"], check=True)
                    index_path = local.with_suffix(".gz.tbi")
                    if not local.is_file() or not index_path.is_file():
                        raise RuntimeError(f"Input staging error for {identifier} {track}.")
                    contigs = set(subprocess.check_output(["tabix", "-l", local], text=True).splitlines())
                    for chrom in chromosomes:
                        present = chrom in contigs
                        if track == "combined" and present:
                            covered.add(chrom)
                        output = slice_root / chrom / "tracks" / f"{local_stem}.{track}.bed.gz"
                        query = ["tabix", "-h", local, chrom] if present else ["tabix", "-H", local]
                        with output.open("wb") as compressed:
                            with subprocess.Popen(query, stdout=subprocess.PIPE) as extracted:
                                try:
                                    subprocess.run(["bgzip", "-c"], stdin=extracted.stdout, stdout=compressed, check=True)
                                finally:
                                    extracted.stdout.close()
                                if extracted.wait() != 0:
                                    raise RuntimeError(f"Tabix failed for {identifier} {track} on {chrom}.")
                        subprocess.run(["tabix", "-p", "bed", output], check=True)
                    local.unlink()
                    index_path.unlink()
            return (identifier, f"tracks/{local_stem}", label), covered

        print(f"Staging {len(samples)} samples for {', '.join(chromosomes)}.", flush=True)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(stage_sample, enumerate(samples, start=1)))

        for chrom in chromosomes:
            coverage_count = sum(chrom in covered for _, covered in results)
            if coverage_count == 0:
                raise SystemExit(f"No combined CpG rows for {chrom} in this subset. Check the contig name and input coverage.")
            chromosome_dir = slice_root / chrom
            with (chromosome_dir / "collection.tsv").open("w", encoding="utf-8") as handle:
                handle.write("identifier\tfilename\tlabels\n")
                for row, _ in results:
                    handle.write("\t".join(row) + "\n")
            with tarfile.open(bundles / f"{chrom}.tar", "w") as archive:
                archive.add(chromosome_dir / "collection.tsv", arcname="collection.tsv")
                archive.add(chromosome_dir / "tracks", arcname="tracks")
            print(f"Staged {chrom}. Samples with combined CpGs: {coverage_count}/{len(samples)}.", flush=True)
        PY
    >>>

    runtime {
        docker: docker_image
        cpu: transfer_workers
        memory: "8 GiB"
        disks: "local-disk ~{disk_gb} SSD"
        maxRetries: 1
    }

    output {
        Array[File]+ chromosome_bundles = glob("bundles/*.tar")
    }
}

task SegmentMethBatChromosome {
    input {
        File chromosome_bundle
        String subset_name
        String chromosome
        Int sample_limit = 0
        Int threads
        Int mem_gb = 32
        Int disk_gb = 50
        Boolean condense_bed_labels
        Boolean enable_nodata_segments
        String docker_image = "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v2"
    }

    String size_suffix = if sample_limit > 0 then ".n~{sample_limit}" else ""
    String result_prefix = "output/" + subset_name + "." + chromosome + size_suffix

    command <<<
        set -euo pipefail
        bundle='~{sub("~{chromosome_bundle}", "'", "'\"'\"'")}'
        subset='~{sub(subset_name, "'", "'\"'\"'")}'
        chrom='~{sub(chromosome, "'", "'\"'\"'")}'
        output_prefix='~{sub(result_prefix, "'", "'\"'\"'")}'
        sample_limit=~{sample_limit}
        if [[ "$bundle" == *://* || ! -r "$bundle" ]]; then
            echo "Segmentation localization error. Cannot read the local chromosome bundle." >&2
            exit 1
        fi
        if [[ ! "$subset" =~ ^subset-[0-9]+$ || ! "$chrom" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
            echo "Invalid subset or chromosome name." >&2
            exit 1
        fi
        if [[ ~{threads} -lt 1 ]]; then
            echo "threads must be at least 1." >&2
            exit 1
        fi
        if [[ "$sample_limit" -lt 0 ]]; then
            echo "sample_limit must be at least 0." >&2
            exit 1
        fi

        echo "Extracting the CpG bundle for $subset on $chrom."
        tar -xf "$bundle"
        if [[ ! -s collection.tsv || ! -d tracks ]]; then
            echo "Chromosome bundle is missing its collection TSV or CpG tracks." >&2
            exit 1
        fi
        awk -v limit="$sample_limit" 'NR == 1 || !limit || NR <= limit + 1' collection.tsv > selected.tsv
        awk -F '\t' 'NR > 1 { print $2 }' selected.tsv > track_prefixes.txt
        sample_count=$(awk 'END {print NR}' track_prefixes.txt)
        if [[ "$sample_count" -lt 1 || ( "$sample_limit" -gt 0 && "$sample_count" -ne "$sample_limit" ) ]]; then
            echo "The bundle does not contain enough samples for $subset on $chrom." >&2
            exit 1
        fi
        while IFS= read -r prefix; do
            for track in combined hap1 hap2; do
                if [[ ! -s "$prefix.$track.bed.gz" || ! -r "$prefix.$track.bed.gz" ]]; then
                    echo "Missing readable CpG track $prefix.$track.bed.gz." >&2
                    exit 1
                fi
            done
        done < track_prefixes.txt

        extra_args=()
        if [[ "~{condense_bed_labels}" == "true" ]]; then
            extra_args+=(--condense-bed-labels)
        fi
        if [[ "~{enable_nodata_segments}" == "true" ]]; then
            extra_args+=(--enable-nodata-segments)
        fi

        mkdir -p output
        echo "Running MethBat joint segmentation for $subset on $chrom with $sample_count samples."
        methbat joint-segment -t ~{threads} --input-collection selected.tsv \
            --output-prefix "$output_prefix" "${extra_args[@]}"
        printf '%s\n' "$sample_count" > sample_count.txt
        printf '%s\t%s\t%s\t%s\n' "$subset" "$chrom" "$sample_limit" "$sample_count" > run_key.tsv
        echo "Finished MethBat joint segmentation for $subset on $chrom."
    >>>

    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{mem_gb} GiB"
        disks: "local-disk ~{disk_gb} SSD"
        maxRetries: 1
    }

    output {
        String result_chromosome = chromosome
        Int sample_count = read_int("sample_count.txt")
        String run_key = read_string("run_key.tsv")
        File meth_regions_bed = result_prefix + ".meth_regions.bed"
        File combined_methyl_bedgraph = result_prefix + ".combined_methyl.bedgraph"
        File asm_bedgraph = result_prefix + ".asm.bedgraph"
    }
}

task CompareMethBatBoundaries {
    input {
        Array[File] region_beds
        Array[String] run_keys
        Float min_reciprocal_overlap
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
    }

    command <<<
        set -euo pipefail

        # Keep each generated list as a File for Cromwell localization.
        bed_list="~{write_lines(prefix("", region_beds))}"
        run_list="~{write_lines(run_keys)}"
        if [[ "$bed_list" == *://* || ! -r "$bed_list" || "$run_list" == *://* || ! -r "$run_list" ]]; then
            echo "Boundary comparison localization error. Cannot read the generated input lists." >&2
            exit 1
        fi
        python3 - --beds "$bed_list" --runs "$run_list" --min-overlap ~{min_reciprocal_overlap} <<'PY'
        import argparse
        from collections import defaultdict
        import csv
        from itertools import combinations
        from math import ceil
        from pathlib import Path
        from statistics import median

        parser = argparse.ArgumentParser()
        parser.add_argument("--beds", type=Path, required=True)
        parser.add_argument("--runs", type=Path, required=True)
        parser.add_argument("--min-overlap", type=float, required=True)
        args = parser.parse_args()

        if not 0 < args.min_overlap <= 1:
            raise SystemExit("min_reciprocal_overlap must be greater than 0 and at most 1.")

        def read_intervals(path, chromosome):
            if "://" in path or not Path(path).is_file():
                raise SystemExit(f"Boundary comparison localization error. Cannot read BED {path}.")
            by_label = defaultdict(list)
            with Path(path).open(encoding="utf-8") as handle:
                for line_number, line in enumerate(handle, start=1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    columns = line.rstrip("\n").split("\t")
                    try:
                        chrom, start, end, label = columns[:4]
                        start, end = int(start), int(end)
                        if chrom != chromosome or start < 0 or end <= start or not label:
                            raise ValueError
                    except ValueError:
                        raise SystemExit(f"Invalid {chromosome} BED interval in {path} line {line_number}.")
                    by_label[label].append((start, end))
            return {label: sorted(intervals) for label, intervals in by_label.items()}

        def reciprocal_best_matches(left, right):
            best_left, best_right = {}, {}
            first_right = 0
            for index_a, (start_a, end_a) in enumerate(left):
                while first_right < len(right) and right[first_right][1] <= start_a:
                    first_right += 1
                for index_b in range(first_right, len(right)):
                    start_b, end_b = right[index_b]
                    if start_b >= end_a:
                        break
                    overlap = min(end_a, end_b) - max(start_a, start_b)
                    if overlap < args.min_overlap * max(end_a - start_a, end_b - start_b):
                        continue
                    union = max(end_a, end_b) - min(start_a, start_b)
                    score = overlap / union
                    # Keep the earlier coordinate on ties.
                    if index_a not in best_left or score > best_left[index_a][0]:
                        best_left[index_a] = (score, index_b)
                    if index_b not in best_right or score > best_right[index_b][0]:
                        best_right[index_b] = (score, index_a)
            return [(i, j) for i, (_, j) in best_left.items() if best_right[j][1] == i]

        def shift_statistics(shifts):
            if not shifts:
                return "NA", "NA"
            values = sorted(abs(shift) for shift in shifts)
            return f"{median(values):g}", values[ceil(0.9 * len(values)) - 1]

        bed_paths = args.beds.read_text(encoding="utf-8").splitlines()
        keys = args.runs.read_text(encoding="utf-8").splitlines()
        if len(bed_paths) != len(keys):
            raise SystemExit("The localized BED list and run metadata have different lengths.")
        by_run = defaultdict(dict)
        for path, key in zip(bed_paths, keys):
            subset, chromosome, sample_limit, sample_count = key.split("\t")
            sample_limit, sample_count = int(sample_limit), int(sample_count)
            runs = by_run[(chromosome, sample_limit)]
            if subset in runs:
                raise SystemExit(f"Duplicate segment output for {subset} on {chromosome} at size {sample_limit}.")
            runs[subset] = (sample_count, read_intervals(path, chromosome))

        fields = ["chromosome", "subset_size", "subset_a", "subset_b", "samples_a", "samples_b", "label"]
        with Path("boundary_summary.tsv").open("w", encoding="utf-8", newline="") as summary_file, \
                Path("boundary_shifts.tsv").open("w", encoding="utf-8", newline="") as shifts_file:
            summary = csv.writer(summary_file, delimiter="\t", lineterminator="\n")
            details = csv.writer(shifts_file, delimiter="\t", lineterminator="\n")
            summary.writerow(fields + [
                "regions_a", "regions_b", "matched_regions", "unmatched_a", "unmatched_b",
                "median_abs_start_shift_bp", "p90_abs_start_shift_bp",
                "median_abs_end_shift_bp", "p90_abs_end_shift_bp",
            ])
            details.writerow(fields + ["start_a", "end_a", "start_b", "end_b", "start_shift_bp", "end_shift_bp"])
            for (chromosome, size), runs in sorted(by_run.items()):
                for subset_a, subset_b in combinations(sorted(runs), 2):
                    count_a, regions_a = runs[subset_a]
                    count_b, regions_b = runs[subset_b]
                    key = [chromosome, size or "all", subset_a, subset_b, count_a, count_b]
                    labels = sorted(regions_a.keys() | regions_b.keys())
                    # Keep a summary row when both runs return an empty BED.
                    for label in labels or [""]:
                        left, right = regions_a.get(label, []), regions_b.get(label, [])
                        matches = reciprocal_best_matches(left, right)
                        starts, ends = [], []
                        for i, j in matches:
                            start_a, end_a = left[i]
                            start_b, end_b = right[j]
                            starts.append(start_b - start_a)
                            ends.append(end_b - end_a)
                            details.writerow(key + [
                                label, start_a, end_a, start_b, end_b, starts[-1], ends[-1],
                            ])
                        n = len(matches)
                        summary.writerow(key + [
                            label, len(left), len(right), n, len(left) - n, len(right) - n,
                            *shift_statistics(starts), *shift_statistics(ends),
                        ])
                    print(f"Compared {subset_a} and {subset_b} on {chromosome}, size {size or 'all'}.", flush=True)
        PY
    >>>

    runtime {
        docker: docker_image
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk 20 HDD"
    }

    output {
        File boundary_summary = "boundary_summary.tsv"
        File boundary_shifts = "boundary_shifts.tsv"
    }
}
