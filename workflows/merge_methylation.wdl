version 1.0

# The global merge is deliberately outside the scatter.  A per-shard site
# filter would use a different denominator in every shard and would therefore
# not implement MinSampleFraction across the cohort.

task ShardMethylationManifest {
    input {
        File SampleManifest
        Int SamplesPerShard
    }

    command <<<
        set -euo pipefail

        if [ "~{SamplesPerShard}" -lt 1 ]; then
            echo "SamplesPerShard must be at least 1" >&2
            exit 1
        fi

        n_samples=$(awk 'END { print NR - 1 }' "~{SampleManifest}")
        if [ "$n_samples" -lt 1 ]; then
            echo "SampleManifest must have a header and at least one sample" >&2
            exit 1
        fi
        printf '%s\n' "$n_samples" > total_samples.txt

        awk -v shard_size="~{SamplesPerShard}" '
            NR == 1 {
                header = $0
                next
            }
            {
                shard = int((NR - 2) / shard_size)
                output = sprintf("methylation_manifest.shard.%05d.tsv", shard)
                if (!(output in seen)) {
                    print header > output
                    seen[output] = 1
                }
                print >> output
            }
        ' "~{SampleManifest}"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "2G"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        Array[File] ShardManifests = glob("methylation_manifest.shard.*.tsv")
        Int TotalSamples = read_int("total_samples.txt")
    }
}

task FilterMethylationShard {
    input {
        File ShardManifest
        String OutputPrefix
        Float MinCoverage
        String FilterChroms
        Float FenceK
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        Rscript /tmp/MergeMethylationCalls.R \
            --InputManifest "~{ShardManifest}" \
            --OutputPrefix "~{OutputPrefix}" \
            --MinCoverage ~{MinCoverage} \
            --FilterChroms "~{FilterChroms}" \
            --FenceK ~{FenceK} \
            --PerSampleOnly
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.per_sample_filtered.long.tsv.gz"
        File AllCalls = "~{OutputPrefix}.methylation.per_sample_qc.long.tsv.gz"
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
    }
}

task MergeMethylationShards {
    input {
        Array[File] FilteredCallShards
        Array[File] AllCallShards
        Array[File] SampleQcShards
        Int TotalSamples
        String OutputPrefix
        Float MinSampleFraction
        Int MinSamples
        String ValueColumn
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' FilteredCallShards} > filtered_call_shards.list
        printf '%s\n' ~{sep=' ' AllCallShards} > all_call_shards.list
        printf '%s\n' ~{sep=' ' SampleQcShards} > sample_qc_shards.list

        Rscript /tmp/MergeMethylationCalls.R \
            --FilteredCallList filtered_call_shards.list \
            --AllCallList all_call_shards.list \
            --FilteredSampleQcList sample_qc_shards.list \
            --TotalSamples ~{TotalSamples} \
            --OutputPrefix "~{OutputPrefix}" \
            --MinSampleFraction ~{MinSampleFraction} \
            --MinSamples ~{MinSamples} \
            --ValueColumn "~{ValueColumn}"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.filtered.long.tsv.gz"
        File SiteQC = "~{OutputPrefix}.methylation.site_qc.tsv.gz"
        File SiteMetadata = "~{OutputPrefix}.methylation.site_metadata.tsv.gz"
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
        File MethylationMatrix = "~{OutputPrefix}.methylation.matrix.bed.gz"
    }
}

workflow MergeMethylation {
    input {
        # TSV with sample_id and absolute file_path values. The paths must be
        # readable by each task container (for example, from a mounted shared
        # filesystem); embedded file references are not localized by WDL.
        File SampleManifest
        String OutputPrefix

        Int SamplesPerShard = 25
        Float MinCoverage = 10.0
        Float MinSampleFraction = 0.8
        Int MinSamples = 0
        String FilterChroms = "X|Y|M|_"
        Float FenceK = 3.0
        String ValueColumn

        Int ShardMemoryGB = 16
        Int ShardDiskGB = 100
        Int MergeMemoryGB = 64
        Int MergeDiskGB = 200
        Int NumThreads = 1
    }

    call ShardMethylationManifest {
        input:
            SampleManifest = SampleManifest,
            SamplesPerShard = SamplesPerShard
    }

    scatter (shard_index in range(length(ShardMethylationManifest.ShardManifests))) {
        File shard_manifest = ShardMethylationManifest.ShardManifests[shard_index]
        String shard_output_prefix = "~{OutputPrefix}.shard.~{shard_index}"

        call FilterMethylationShard {
            input:
                ShardManifest = shard_manifest,
                OutputPrefix = shard_output_prefix,
                MinCoverage = MinCoverage,
                FilterChroms = FilterChroms,
                FenceK = FenceK,
                MemoryGB = ShardMemoryGB,
                DiskGB = ShardDiskGB,
                NumThreads = NumThreads
        }
    }

    call MergeMethylationShards {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCalls,
            AllCallShards = FilterMethylationShard.AllCalls,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            OutputPrefix = OutputPrefix,
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            ValueColumn = ValueColumn,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    output {
        File FilteredCalls = MergeMethylationShards.FilteredCalls
        File SiteQC = MergeMethylationShards.SiteQC
        File SiteMetadata = MergeMethylationShards.SiteMetadata
        File SampleQC = MergeMethylationShards.SampleQC
        File MethylationMatrix = MergeMethylationShards.MethylationMatrix
        Array[File] ShardSampleQC = FilterMethylationShard.SampleQC
    }
}
