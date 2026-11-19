version 1.0

workflow MethBat_jointseg_chromosome_subsets {
    meta {
        description: "Split a pb-CpG-tools cohort into disjoint subsets and run MethBat v0.17 joint segmentation by chromosome."
    }

    input {
        File manifest_file
        Array[String]+ chromosomes = ["chr21"]
        Int subset_count = 2
        Int split_seed = 42

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
        subset_count: "Number of non-overlapping, label-balanced subsets. Use 1 for the full cohort."
        split_seed: "Random seed for reproducible subset membership."
    }

    call SplitMethBatCohort {
        input:
            manifest_file = manifest_file,
            subset_count = subset_count,
            split_seed = split_seed,
            docker_image = staging_docker
    }

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

        scatter (bundle in StageMethBatChromosomes.chromosome_bundles) {
            call SegmentMethBatChromosome {
                input:
                    chromosome_bundle = bundle,
                    subset_name = subset_name,
                    chromosome = basename(bundle, ".tar"),
                    threads = threads,
                    mem_gb = mem_gb,
                    disk_gb = disk_gb,
                    condense_bed_labels = condense_bed_labels,
                    enable_nodata_segments = enable_nodata_segments,
                    docker_image = methbat_docker
            }
        }
    }

    output {
        File membership_tsv = SplitMethBatCohort.membership_tsv
        Array[File] subset_manifests = SplitMethBatCohort.subset_manifests
        # Arrays are indexed first by subset and then by chromosome.
        Array[String] subsets = subset_name
        Array[Array[String]] result_chromosomes = SegmentMethBatChromosome.result_chromosome
        Array[Array[File]] meth_regions_bed = SegmentMethBatChromosome.meth_regions_bed
        Array[Array[File]] combined_methyl_bedgraph = SegmentMethBatChromosome.combined_methyl_bedgraph
        Array[Array[File]] asm_bedgraph = SegmentMethBatChromosome.asm_bedgraph
    }
}

task SplitMethBatCohort {
    input {
        File manifest_file
        Int subset_count
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

        python3 - --manifest "$manifest" --subsets ~{subset_count} --seed ~{split_seed} <<'PY'
        import argparse
        from collections import defaultdict
        from pathlib import Path
        import random
        import re

        parser = argparse.ArgumentParser()
        parser.add_argument("--manifest", type=Path, required=True)
        parser.add_argument("--subsets", type=int, required=True)
        parser.add_argument("--seed", type=int, required=True)
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

        groups = [[] for _ in range(args.subsets)]
        rng = random.Random(args.seed)
        offset = 0
        for label in sorted(samples_by_label):
            samples = sorted(samples_by_label[label])
            rng.shuffle(samples)
            for index, sample in enumerate(samples):
                groups[(index + offset) % args.subsets].append(sample)
            offset = (offset + len(samples)) % args.subsets

        output = Path("subsets")
        output.mkdir()
        with Path("membership.tsv").open("w", encoding="utf-8") as membership:
            membership.write("identifier\tlabels\tsubset\n")
            for number, group in enumerate(groups, start=1):
                subset = f"subset-{number:06d}"
                with (output / f"{subset}.tsv").open("w", encoding="utf-8") as handle:
                    handle.write("identifier\tfilename\tlabels\n")
                    for identifier, prefix, label in sorted(group):
                        handle.write(f"{identifier}\t{prefix}\t{label}\n")
                        membership.write(f"{identifier}\t{label}\t{subset}\n")
                print(f"Wrote {subset} with {len(group)} samples.", flush=True)
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
        chromosome_list='~{sub("~{write_lines(chromosomes)}", "'", "'\"'\"'")}'
        if [[ "$manifest" == *://* || ! -r "$manifest" || ! -r "$chromosome_list" ]]; then
            echo "Staging localization error. Cannot read the subset manifest or chromosome list." >&2
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
        Int threads
        Int mem_gb = 32
        Int disk_gb = 50
        Boolean condense_bed_labels
        Boolean enable_nodata_segments
        String docker_image = "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v2"
    }

    String result_prefix = "output/" + subset_name + "." + chromosome

    command <<<
        set -euo pipefail
        bundle='~{sub("~{chromosome_bundle}", "'", "'\"'\"'")}'
        subset='~{sub(subset_name, "'", "'\"'\"'")}'
        chrom='~{sub(chromosome, "'", "'\"'\"'")}'
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

        echo "Extracting the CpG bundle for $subset on $chrom."
        tar -xf "$bundle"
        if [[ ! -s collection.tsv || ! -d tracks ]]; then
            echo "Chromosome bundle is missing its collection TSV or CpG tracks." >&2
            exit 1
        fi
        awk -F '\t' 'NR > 1 { print $2 }' collection.tsv > track_prefixes.txt
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
        echo "Running MethBat joint segmentation for $subset on $chrom."
        methbat joint-segment -t ~{threads} --input-collection collection.tsv \
            --output-prefix "output/$subset.$chrom" "${extra_args[@]}"
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
        File meth_regions_bed = result_prefix + ".meth_regions.bed"
        File combined_methyl_bedgraph = result_prefix + ".combined_methyl.bedgraph"
        File asm_bedgraph = result_prefix + ".asm.bedgraph"
    }
}
