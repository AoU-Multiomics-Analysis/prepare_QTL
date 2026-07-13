version 1.0
import "calculate_phenotypePCs.wdl" as ComputePCs
import "MergeCovariates.wdl" as CovariateMerge

# Cohort-level entry point for per-sample methylation outputs. Terra receives
# one compact manifest File rather than 22 large Array[File] inputs.

task PrepareMethylationCohortManifest {
    input {
        File CohortManifest
    }

    command <<<
        Rscript /tmp/PrepareMethylationCohortManifest.R \
            --CohortManifest "~{CohortManifest}" \
            --OutputDir manifest_lists
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "4G"
        disks: "local-disk 20 HDD"
        cpu: 1
    }

    output {
        File SampleQCManifest = "manifest_lists/sample_qc_paths.list"
        File AllCallsAutosome01Manifest = "manifest_lists/autosome01_paths.list"
        File AllCallsAutosome02Manifest = "manifest_lists/autosome02_paths.list"
        File AllCallsAutosome03Manifest = "manifest_lists/autosome03_paths.list"
        File AllCallsAutosome04Manifest = "manifest_lists/autosome04_paths.list"
        File AllCallsAutosome05Manifest = "manifest_lists/autosome05_paths.list"
        File AllCallsAutosome06Manifest = "manifest_lists/autosome06_paths.list"
        File AllCallsAutosome07Manifest = "manifest_lists/autosome07_paths.list"
        File AllCallsAutosome08Manifest = "manifest_lists/autosome08_paths.list"
        File AllCallsAutosome09Manifest = "manifest_lists/autosome09_paths.list"
        File AllCallsAutosome10Manifest = "manifest_lists/autosome10_paths.list"
        File AllCallsAutosome11Manifest = "manifest_lists/autosome11_paths.list"
        File AllCallsAutosome12Manifest = "manifest_lists/autosome12_paths.list"
        File AllCallsAutosome13Manifest = "manifest_lists/autosome13_paths.list"
        File AllCallsAutosome14Manifest = "manifest_lists/autosome14_paths.list"
        File AllCallsAutosome15Manifest = "manifest_lists/autosome15_paths.list"
        File AllCallsAutosome16Manifest = "manifest_lists/autosome16_paths.list"
        File AllCallsAutosome17Manifest = "manifest_lists/autosome17_paths.list"
        File AllCallsAutosome18Manifest = "manifest_lists/autosome18_paths.list"
        File AllCallsAutosome19Manifest = "manifest_lists/autosome19_paths.list"
        File AllCallsAutosome20Manifest = "manifest_lists/autosome20_paths.list"
        File AllCallsAutosome21Manifest = "manifest_lists/autosome21_paths.list"
        File AllCallsAutosome22Manifest = "manifest_lists/autosome22_paths.list"
    }
}

task BuildMethylationCohortSamples {
    input {
        Array[File] SampleQCFiles
        String OutputPrefix
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' SampleQCFiles} > sample_qc_files.list
        Rscript /tmp/BuildMethylationCohortSamples.R \
            --SampleQcList sample_qc_files.list \
            --OutputPrefix "~{OutputPrefix}"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "4G"
        disks: "local-disk 20 HDD"
        cpu: 1
    }

    output {
        File CohortSamples = "~{OutputPrefix}.methylation.cohort_samples.tsv"
        File CohortSampleQC = "~{OutputPrefix}.methylation.cohort_sample_qc.tsv"
        Int TotalSamples = read_int("~{OutputPrefix}.methylation.total_samples.txt")
    }
}

task MergeMethylationChromosome {
    input {
        Array[File] AllCallShards
        File CohortSampleQC
        File CohortSamples
        Int TotalSamples
        String Chromosome
        String OutputPrefix
        Float MinSampleFraction
        Int MinSamples
        Float MinMethylationMAD
        String ValueColumn
        Float ValueMultiplier
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' AllCallShards} > all_call_shards.list
        methylation-chromosome-merge \
            --all-call-list all_call_shards.list \
            --sample-qc "~{CohortSampleQC}" \
            --cohort-samples "~{CohortSamples}" \
            --total-samples ~{TotalSamples} \
            --chromosome "~{Chromosome}" \
            --output-prefix "~{OutputPrefix}" \
            --min-sample-fraction ~{MinSampleFraction} \
            --min-samples ~{MinSamples} \
            --min-methylation-mad ~{MinMethylationMAD} \
            --value-column "~{ValueColumn}" \
            --value-multiplier ~{ValueMultiplier} \
            --progress-every-sites 1000
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl-methylation-rust:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.filtered.long.tsv.gz"
        File SiteQC = "~{OutputPrefix}.methylation.site_qc.tsv.gz"
        File SiteMetadata = "~{OutputPrefix}.methylation.site_metadata.tsv.gz"
        File RawMethylationBed = "~{OutputPrefix}.methylation.raw.bed.gz"
        File IntMethylationBed = "~{OutputPrefix}.methylation.INT.bed.gz"
    }
}

task BuildMethylationCorrelationCovariates {
    input {
        File PhenotypePCs
        File? AdditionalCovariates
        String OutputPrefix
    }

    command <<<
        set -euo pipefail
        Rscript /tmp/BuildMethylationCorrelationCovariates.R \
            --PhenotypePCs "~{PhenotypePCs}" \
            ~{if defined(AdditionalCovariates) then "--AdditionalCovariates '" + select_first([AdditionalCovariates]) + "'" else ""} \
            --OutputFile "~{OutputPrefix}.methylation.correlation_covariates.tsv"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "16G"
        disks: "local-disk 50 HDD"
        cpu: 1
    }

    output {
        File CorrelationCovariates = "~{OutputPrefix}.methylation.correlation_covariates.tsv"
    }
}

task AnalyzeMethylationCpGCorrelation {
    input {
        File IntMethylationBed
        File Covariates
        String OutputPrefix
        Int WindowBP
        Float MinAbsCorrelation
        Int MemoryGB
        Int DiskGB
    }

    command <<<
        Rscript /tmp/AnalyzeMethylationCpGCorrelation.R \
            --InputBed "~{IntMethylationBed}" \
            --Covariates "~{Covariates}" \
            --OutputPrefix "~{OutputPrefix}" \
            --WindowBP ~{WindowBP} \
            --MinAbsCorrelation ~{MinAbsCorrelation}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File CorrelationClusters = "~{OutputPrefix}.methylation.cpg_correlation_clusters.tsv.gz"
        File RepresentativeCpGs = "~{OutputPrefix}.methylation.cpg_correlation_representatives.tsv.gz"
        File CorrelationSummary = "~{OutputPrefix}.methylation.cpg_correlation_summary.tsv"
        File ClusterSizePlot = "~{OutputPrefix}.methylation.cpg_correlation_cluster_sizes.png"
        File MeanDistancePlot = "~{OutputPrefix}.methylation.cpg_correlation_mean_distances.png"
        File MaxDistancePlot = "~{OutputPrefix}.methylation.cpg_correlation_max_distances.png"
        File SpanPlot = "~{OutputPrefix}.methylation.cpg_correlation_span_vs_size.png"
    }
}

task AggregateMethylationChromosomes {
    input {
        Array[File] FilteredCallsByChromosome
        Array[File] SiteQCByChromosome
        Array[File] SiteMetadataByChromosome
        Array[File] RawMethylationBedByChromosome
        Array[File] IntMethylationBedByChromosome
        File CohortSampleQC
        Int TotalSamples
        String OutputPrefix
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' FilteredCallsByChromosome} > filtered_calls_by_chromosome.list
        printf '%s\n' ~{sep=' ' SiteQCByChromosome} > site_qc_by_chromosome.list
        printf '%s\n' ~{sep=' ' SiteMetadataByChromosome} > site_metadata_by_chromosome.list
        printf '%s\n' ~{sep=' ' RawMethylationBedByChromosome} > raw_beds_by_chromosome.list
        printf '%s\n' ~{sep=' ' IntMethylationBedByChromosome} > int_beds_by_chromosome.list
        printf '%s\n' "~{CohortSampleQC}" > sample_qc_files.list

        concat_chromosome_tables() {
            local list_path="$1"
            local output_path="$2"
            local label="$3"
            local expected_header=""
            local input_path
            local current_header

            while IFS= read -r input_path; do
                [ -n "$input_path" ] || continue
                current_header=$(zgrep -m 1 '^' "$input_path")
                if [ -z "$expected_header" ]; then
                    expected_header="$current_header"
                elif [ "$current_header" != "$expected_header" ]; then
                    echo "${label} chromosome files do not have identical headers: ${input_path}" >&2
                    exit 1
                fi
            done < "$list_path"

            if [ -z "$expected_header" ]; then
                echo "${label} chromosome file list is empty" >&2
                exit 1
            fi

            {
                local first_file=1
                while IFS= read -r input_path; do
                    [ -n "$input_path" ] || continue
                    if [ "$first_file" -eq 1 ]; then
                        bgzip -c -d -@ ~{NumThreads} "$input_path"
                        first_file=0
                    else
                        bgzip -c -d -@ ~{NumThreads} "$input_path" | tail -n +2
                    fi
                done < "$list_path"
            } | bgzip -c -@ ~{NumThreads} > "$output_path"
        }

        concat_chromosome_tables filtered_calls_by_chromosome.list \
            "~{OutputPrefix}.methylation.filtered.long.pre_connectivity.tsv.gz" "Filtered-call"
        concat_chromosome_tables site_qc_by_chromosome.list \
            "~{OutputPrefix}.methylation.site_qc.tsv.gz" "Site-QC"
        concat_chromosome_tables site_metadata_by_chromosome.list \
            "~{OutputPrefix}.methylation.site_metadata.tsv.gz" "Site-metadata"
        concat_chromosome_tables raw_beds_by_chromosome.list \
            "~{OutputPrefix}.methylation.raw.pre_connectivity.bed.gz" "Raw BED"
        concat_chromosome_tables int_beds_by_chromosome.list \
            "~{OutputPrefix}.methylation.INT.pre_connectivity.bed.gz" "INT BED"

        extract_passing_site_metadata() {
            local input_path="$1"
            local output_path="$2"

            bgzip -c -d -@ ~{NumThreads} "$input_path" | \
                awk 'BEGIN { FS = OFS = "\t" }
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            column_name = $i
                            gsub(/^"|"$/, "", column_name)
                            if (column_name == "keep_site") keep_column = i
                        }
                        if (!keep_column) {
                            print "Site-metadata table is missing keep_site" > "/dev/stderr"
                            exit 1
                        }
                        print
                        next
                    }
                    {
                        keep_value = $keep_column
                        gsub(/^"|"$/, "", keep_value)
                        if (keep_value == "TRUE") print
                    }' | \
                bgzip -c -@ ~{NumThreads} > "$output_path"
        }

        extract_passing_site_metadata \
            "~{OutputPrefix}.methylation.site_metadata.tsv.gz" \
            "~{OutputPrefix}.methylation.passing_site_metadata.tsv.gz"

        Rscript /tmp/AggregateMethylationChromosomes.R \
            --SiteMetadata "~{OutputPrefix}.methylation.site_metadata.tsv.gz" \
            --SampleQcList sample_qc_files.list \
            --TotalSamples ~{TotalSamples} \
            --OutputPrefix "~{OutputPrefix}" \
            --SampleQcOutput "~{OutputPrefix}.methylation.sample_qc.pre_connectivity.tsv"

    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File PreConnectivityFilteredCalls = "~{OutputPrefix}.methylation.filtered.long.pre_connectivity.tsv.gz"
        File SiteQC = "~{OutputPrefix}.methylation.site_qc.tsv.gz"
        File SiteMetadata = "~{OutputPrefix}.methylation.site_metadata.tsv.gz"
        File PassingSiteMetadata = "~{OutputPrefix}.methylation.passing_site_metadata.tsv.gz"
        File PreConnectivitySampleQC = "~{OutputPrefix}.methylation.sample_qc.pre_connectivity.tsv"
        File FilterSummary = "~{OutputPrefix}.methylation.filter_summary.tsv"
        File FilterCountsPlot = "~{OutputPrefix}.methylation.filter_counts.png"
        File FilterUpsetPlot = "~{OutputPrefix}.methylation.filter_upset.png"
        File PreConnectivityRawMethylationBed = "~{OutputPrefix}.methylation.raw.pre_connectivity.bed.gz"
        File PreConnectivityIntMethylationBed = "~{OutputPrefix}.methylation.INT.pre_connectivity.bed.gz"
    }
}

task FinalizeMethylationConnectivity {
    input {
        Array[File] RepresentativeCpGsByChromosome
        File PreConnectivityFilteredCalls
        File PreConnectivityRawMethylationBed
        File PreConnectivityIntMethylationBed
        File PreConnectivitySampleQC
        String OutputPrefix
        Int MaxConnectivityFeatures
        Int ConnectivityLandmarks
        Float ConnectivityZThreshold
        Int MemoryGB
        Int DiskGB
    }

    command <<<
        set -euo pipefail
        printf '%s\n' "~{PreConnectivityIntMethylationBed}" > int_beds_by_chromosome.list
        printf '%s\n' ~{sep=' ' RepresentativeCpGsByChromosome} > representative_cpgs_by_chromosome.list

        Rscript /tmp/FinalizeMethylationConnectivity.R \
            --IntBedList int_beds_by_chromosome.list \
            --RepresentativeList representative_cpgs_by_chromosome.list \
            --FilteredCalls "~{PreConnectivityFilteredCalls}" \
            --RawBed "~{PreConnectivityRawMethylationBed}" \
            --IntBed "~{PreConnectivityIntMethylationBed}" \
            --SampleQC "~{PreConnectivitySampleQC}" \
            --OutputPrefix "~{OutputPrefix}" \
            --MaxConnectivityFeatures ~{MaxConnectivityFeatures} \
            --ConnectivityLandmarks ~{ConnectivityLandmarks} \
            --ConnectivityZThreshold ~{ConnectivityZThreshold}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.filtered.long.tsv.gz"
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
        File RawMethylationBed = "~{OutputPrefix}.methylation.raw.bed.gz"
        File IntMethylationBed = "~{OutputPrefix}.methylation.INT.bed.gz"
        File ConnectivityOutliers = "~{OutputPrefix}.methylation.connectivity_outliers.tsv"
        File ConnectivitySummary = "~{OutputPrefix}.methylation.connectivity_summary.tsv"
        File ConnectivityRepresentativeCpGs = "~{OutputPrefix}.methylation.connectivity_representative_cpgs.tsv.gz"
    }
}

task AnnotateMethylationSites {
    input {
        File PassingSiteMetadata
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations
        String OutputPrefix
        Int PromoterWindow
        Int MemoryGB
        Int DiskGB
    }

    command <<<
        Rscript /tmp/AnnotateMethylationSites.R \
            --PassingSiteMetadata "~{PassingSiteMetadata}" \
            --AnnotationGTF "~{AnnotationGTF}" \
            --CCREAnnotations "~{CCREAnnotations}" \
            --CpGIslandAnnotations "~{CpGIslandAnnotations}" \
            --OutputPrefix "~{OutputPrefix}" \
            --PromoterWindow ~{PromoterWindow}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File PassingSiteAnnotations = "~{OutputPrefix}.methylation.passing_site_annotations.tsv.gz"
    }
}

workflow AggregateMethylationCohort {
    input {
        File CohortManifest
        String OutputPrefix
        File? AdditionalCovariates
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations

        Float MinSampleFraction = 0.95
        Int MinSamples = 0
        Float MinMethylationMAD = 0.003
        String AutosomePrefix = "chr"
        Int PromoterWindow = 2000
        String ValueColumn = "mod_score"
        Float ValueMultiplier = 0.01
        Int MergeMemoryGB = 128
        Int MergeDiskGB = 500
        Int AggregateMemoryGB = 64
        Int AggregateDiskGB = 1000
        Int AnnotationMemoryGB = 64
        Int AnnotationDiskGB = 100
        Int CorrelationWindowBP = 1000
        Float CorrelationMinAbsCorrelation = 0.95
        Int CorrelationMemoryGB = 64
        Int CorrelationDiskGB = 250
        Int MaxConnectivityFeatures = 0
        Int ConnectivityLandmarks = 200
        Float ConnectivityZThreshold = -3.0
        Int NumThreads = 1
    }

    call PrepareMethylationCohortManifest {
        input:
            CohortManifest = CohortManifest
    }

    Array[File] SampleQCFiles = read_lines(PrepareMethylationCohortManifest.SampleQCManifest)
    Array[Array[File]] AllCallFilesByAutosome = [
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome01Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome02Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome03Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome04Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome05Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome06Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome07Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome08Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome09Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome10Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome11Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome12Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome13Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome14Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome15Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome16Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome17Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome18Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome19Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome20Manifest),
        read_lines(PrepareMethylationCohortManifest.AllCallsAutosome21Manifest), read_lines(PrepareMethylationCohortManifest.AllCallsAutosome22Manifest)
    ]

    call BuildMethylationCohortSamples {
        input:
            SampleQCFiles = SampleQCFiles,
            OutputPrefix = OutputPrefix
    }

    Array[String] AutosomeNames = [
        AutosomePrefix + "1", AutosomePrefix + "2", AutosomePrefix + "3", AutosomePrefix + "4",
        AutosomePrefix + "5", AutosomePrefix + "6", AutosomePrefix + "7", AutosomePrefix + "8",
        AutosomePrefix + "9", AutosomePrefix + "10", AutosomePrefix + "11", AutosomePrefix + "12",
        AutosomePrefix + "13", AutosomePrefix + "14", AutosomePrefix + "15", AutosomePrefix + "16",
        AutosomePrefix + "17", AutosomePrefix + "18", AutosomePrefix + "19", AutosomePrefix + "20",
        AutosomePrefix + "21", AutosomePrefix + "22"
    ]
    Array[String] AutosomeOutputSuffixes = [
        "autosome01", "autosome02", "autosome03", "autosome04",
        "autosome05", "autosome06", "autosome07", "autosome08",
        "autosome09", "autosome10", "autosome11", "autosome12",
        "autosome13", "autosome14", "autosome15", "autosome16",
        "autosome17", "autosome18", "autosome19", "autosome20",
        "autosome21", "autosome22"
    ]
    scatter (autosome_index in range(length(AutosomeNames))) {
        call MergeMethylationChromosome as MergeMethylationAutosome {
            input:
                AllCallShards = AllCallFilesByAutosome[autosome_index],
                CohortSampleQC = BuildMethylationCohortSamples.CohortSampleQC,
                CohortSamples = BuildMethylationCohortSamples.CohortSamples,
                TotalSamples = BuildMethylationCohortSamples.TotalSamples,
                Chromosome = AutosomeNames[autosome_index],
                OutputPrefix = OutputPrefix + "." + AutosomeOutputSuffixes[autosome_index],
                MinSampleFraction = MinSampleFraction,
                MinSamples = MinSamples,
                MinMethylationMAD = MinMethylationMAD,
                ValueColumn = ValueColumn,
                ValueMultiplier = ValueMultiplier,
                MemoryGB = MergeMemoryGB,
                DiskGB = MergeDiskGB,
                NumThreads = NumThreads
        }
    }

    call AggregateMethylationChromosomes {
        input:
            FilteredCallsByChromosome = MergeMethylationAutosome.FilteredCalls,
            SiteQCByChromosome = MergeMethylationAutosome.SiteQC,
            SiteMetadataByChromosome = MergeMethylationAutosome.SiteMetadata,
            RawMethylationBedByChromosome = MergeMethylationAutosome.RawMethylationBed,
            IntMethylationBedByChromosome = MergeMethylationAutosome.IntMethylationBed,
            CohortSampleQC = BuildMethylationCohortSamples.CohortSampleQC,
            TotalSamples = BuildMethylationCohortSamples.TotalSamples,
            OutputPrefix = OutputPrefix,
            MemoryGB = AggregateMemoryGB,
            DiskGB = AggregateDiskGB,
            NumThreads = NumThreads
    }

    call ComputePCs.PhenotypePCs as PreliminaryIntPhenotypePCs {
        input:
            BedFile = AggregateMethylationChromosomes.PreConnectivityIntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation.pre_connectivity",
            OutputSuffix = ".INT",
            memory = MergeMemoryGB,
            disk_space = MergeDiskGB,
            num_threads = NumThreads
    }

    call BuildMethylationCorrelationCovariates {
        input:
            PhenotypePCs = PreliminaryIntPhenotypePCs.OutPhenotypePCs,
            AdditionalCovariates = AdditionalCovariates,
            OutputPrefix = OutputPrefix
    }

    scatter (autosome_index in range(length(AutosomeNames))) {
        call AnalyzeMethylationCpGCorrelation as AnalyzeMethylationAutosomeCorrelation {
            input:
                IntMethylationBed = MergeMethylationAutosome.IntMethylationBed[autosome_index],
                Covariates = BuildMethylationCorrelationCovariates.CorrelationCovariates,
                OutputPrefix = OutputPrefix + "." + AutosomeOutputSuffixes[autosome_index],
                WindowBP = CorrelationWindowBP,
                MinAbsCorrelation = CorrelationMinAbsCorrelation,
                MemoryGB = CorrelationMemoryGB,
                DiskGB = CorrelationDiskGB
        }
    }

    call FinalizeMethylationConnectivity {
        input:
            RepresentativeCpGsByChromosome = AnalyzeMethylationAutosomeCorrelation.RepresentativeCpGs,
            PreConnectivityFilteredCalls = AggregateMethylationChromosomes.PreConnectivityFilteredCalls,
            PreConnectivityRawMethylationBed = AggregateMethylationChromosomes.PreConnectivityRawMethylationBed,
            PreConnectivityIntMethylationBed = AggregateMethylationChromosomes.PreConnectivityIntMethylationBed,
            PreConnectivitySampleQC = AggregateMethylationChromosomes.PreConnectivitySampleQC,
            OutputPrefix = OutputPrefix,
            MaxConnectivityFeatures = MaxConnectivityFeatures,
            ConnectivityLandmarks = ConnectivityLandmarks,
            ConnectivityZThreshold = ConnectivityZThreshold,
            MemoryGB = AggregateMemoryGB,
            DiskGB = AggregateDiskGB
    }

    call AnnotateMethylationSites {
        input:
            PassingSiteMetadata = AggregateMethylationChromosomes.PassingSiteMetadata,
            AnnotationGTF = AnnotationGTF,
            CCREAnnotations = CCREAnnotations,
            CpGIslandAnnotations = CpGIslandAnnotations,
            OutputPrefix = OutputPrefix,
            PromoterWindow = PromoterWindow,
            MemoryGB = AnnotationMemoryGB,
            DiskGB = AnnotationDiskGB
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = FinalizeMethylationConnectivity.IntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".INT",
            memory = MergeMemoryGB,
            disk_space = MergeDiskGB,
            num_threads = NumThreads
    }

    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".INT"
        }
    }

    output {
        File FilteredCalls = FinalizeMethylationConnectivity.FilteredCalls
        File SiteQC = AggregateMethylationChromosomes.SiteQC
        File SiteMetadata = AggregateMethylationChromosomes.SiteMetadata
        File PassingSiteMetadata = AggregateMethylationChromosomes.PassingSiteMetadata
        File SampleQC = FinalizeMethylationConnectivity.SampleQC
        File FilterSummary = AggregateMethylationChromosomes.FilterSummary
        File FilterCountsPlot = AggregateMethylationChromosomes.FilterCountsPlot
        File FilterUpsetPlot = AggregateMethylationChromosomes.FilterUpsetPlot
        File RawMethylationBed = FinalizeMethylationConnectivity.RawMethylationBed
        File IntMethylationBed = FinalizeMethylationConnectivity.IntMethylationBed
        File ConnectivityOutliers = FinalizeMethylationConnectivity.ConnectivityOutliers
        File ConnectivitySummary = FinalizeMethylationConnectivity.ConnectivitySummary
        File ConnectivityRepresentativeCpGs = FinalizeMethylationConnectivity.ConnectivityRepresentativeCpGs
        Array[File] CorrelationClustersByChromosome = AnalyzeMethylationAutosomeCorrelation.CorrelationClusters
        Array[File] CorrelationSummariesByChromosome = AnalyzeMethylationAutosomeCorrelation.CorrelationSummary
        File PassingSiteAnnotations = AnnotateMethylationSites.PassingSiteAnnotations
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File CohortSamples = BuildMethylationCohortSamples.CohortSamples
        Int TotalSamples = BuildMethylationCohortSamples.TotalSamples
    }
}
