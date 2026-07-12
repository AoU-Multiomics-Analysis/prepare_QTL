version 1.0

# Terra-table entry point: run once per sample entity, with SampleID and
# MethylationBed bound directly to columns on that table. No external manifest
# or shared input staging is required.

task FilterMethylationSample {
    input {
        String SampleID
        File MethylationBed
        String OutputPrefix
        Float MinCoverage
        String FilterChroms
        Float FenceK
        String AutosomePrefix
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail

        if [[ ! "~{SampleID}" =~ ^[A-Za-z0-9._-]+$ ]]; then
            echo "SampleID must match [A-Za-z0-9._-]+" >&2
            exit 1
        fi
        if [ ~{NumThreads} -lt 1 ]; then
            echo "NumThreads must be at least 1" >&2
            exit 1
        fi

        # The supplied File is localized by Cromwell. This one-row manifest is
        # task-local implementation detail only; users do not provide one.
        printf 'sample_id\tfile_path\n%s\t%s\n' "~{SampleID}" "~{MethylationBed}" > input_manifest.tsv

        Rscript /opt/prepare_qtl/scripts/methylation/FilterMethylationShard.R \
            --InputManifest input_manifest.tsv \
            --OutputPrefix "~{OutputPrefix}" \
            --MinCoverage ~{MinCoverage} \
            --FilterChroms "~{FilterChroms}" \
            --FenceK ~{FenceK} \
            --AutosomePrefix "~{AutosomePrefix}" \
            --NumThreads ~{NumThreads}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
        File AllCallsAutosome01 = "~{OutputPrefix}.methylation.autosome01.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome02 = "~{OutputPrefix}.methylation.autosome02.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome03 = "~{OutputPrefix}.methylation.autosome03.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome04 = "~{OutputPrefix}.methylation.autosome04.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome05 = "~{OutputPrefix}.methylation.autosome05.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome06 = "~{OutputPrefix}.methylation.autosome06.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome07 = "~{OutputPrefix}.methylation.autosome07.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome08 = "~{OutputPrefix}.methylation.autosome08.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome09 = "~{OutputPrefix}.methylation.autosome09.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome10 = "~{OutputPrefix}.methylation.autosome10.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome11 = "~{OutputPrefix}.methylation.autosome11.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome12 = "~{OutputPrefix}.methylation.autosome12.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome13 = "~{OutputPrefix}.methylation.autosome13.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome14 = "~{OutputPrefix}.methylation.autosome14.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome15 = "~{OutputPrefix}.methylation.autosome15.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome16 = "~{OutputPrefix}.methylation.autosome16.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome17 = "~{OutputPrefix}.methylation.autosome17.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome18 = "~{OutputPrefix}.methylation.autosome18.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome19 = "~{OutputPrefix}.methylation.autosome19.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome20 = "~{OutputPrefix}.methylation.autosome20.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome21 = "~{OutputPrefix}.methylation.autosome21.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome22 = "~{OutputPrefix}.methylation.autosome22.per_sample_qc.long.tsv.gz"
    }
}

workflow ProcessMethylationSample {
    input {
        String SampleID
        File MethylationBed
        Float MinCoverage = 10.0
        String FilterChroms = "X|Y|M|_"
        String AutosomePrefix = "chr"
        Float FenceK = 3.0
        Int MemoryGB = 64
        Int DiskGB = 250
        Int NumThreads = 4
    }

    call FilterMethylationSample {
        input:
            SampleID = SampleID,
            MethylationBed = MethylationBed,
            OutputPrefix = SampleID,
            MinCoverage = MinCoverage,
            FilterChroms = FilterChroms,
            FenceK = FenceK,
            AutosomePrefix = AutosomePrefix,
            MemoryGB = MemoryGB,
            DiskGB = DiskGB,
            NumThreads = NumThreads
    }

    output {
        File SampleQC = FilterMethylationSample.SampleQC
        File AllCallsAutosome01 = FilterMethylationSample.AllCallsAutosome01
        File AllCallsAutosome02 = FilterMethylationSample.AllCallsAutosome02
        File AllCallsAutosome03 = FilterMethylationSample.AllCallsAutosome03
        File AllCallsAutosome04 = FilterMethylationSample.AllCallsAutosome04
        File AllCallsAutosome05 = FilterMethylationSample.AllCallsAutosome05
        File AllCallsAutosome06 = FilterMethylationSample.AllCallsAutosome06
        File AllCallsAutosome07 = FilterMethylationSample.AllCallsAutosome07
        File AllCallsAutosome08 = FilterMethylationSample.AllCallsAutosome08
        File AllCallsAutosome09 = FilterMethylationSample.AllCallsAutosome09
        File AllCallsAutosome10 = FilterMethylationSample.AllCallsAutosome10
        File AllCallsAutosome11 = FilterMethylationSample.AllCallsAutosome11
        File AllCallsAutosome12 = FilterMethylationSample.AllCallsAutosome12
        File AllCallsAutosome13 = FilterMethylationSample.AllCallsAutosome13
        File AllCallsAutosome14 = FilterMethylationSample.AllCallsAutosome14
        File AllCallsAutosome15 = FilterMethylationSample.AllCallsAutosome15
        File AllCallsAutosome16 = FilterMethylationSample.AllCallsAutosome16
        File AllCallsAutosome17 = FilterMethylationSample.AllCallsAutosome17
        File AllCallsAutosome18 = FilterMethylationSample.AllCallsAutosome18
        File AllCallsAutosome19 = FilterMethylationSample.AllCallsAutosome19
        File AllCallsAutosome20 = FilterMethylationSample.AllCallsAutosome20
        File AllCallsAutosome21 = FilterMethylationSample.AllCallsAutosome21
        File AllCallsAutosome22 = FilterMethylationSample.AllCallsAutosome22
    }
}
