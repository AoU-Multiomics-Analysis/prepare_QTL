version 1.0
import "cohort_aggregation.wdl" as CohortAggregation
import "connectivity.wdl" as Connectivity
import "qtl_covariates.wdl" as QtlCovariates

# Public cohort-level entry point for per-sample methylation outputs. The input
# and output contract stays stable while implementation is split into logical
# stage workflows in this directory.

workflow AggregateMethylationCohort {
    input {
        String annotation__annotate_methylation_sites_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String calculate_phenotypepcs__computep_cs_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String cohort_aggregation__aggregate_methylation_chromosomes_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String cohort_aggregation__build_methylation_cohort_samples_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String cohort_aggregation__merge_methylation_chromosome_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-methylation-rust@sha256:16f631c34e0ce265d686335b91c18948607127178d95e7829070c97cd207d6ad"
        String cohort_aggregation__prepare_methylation_cohort_manifest_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String connectivity__analyze_methylation_cpg_correlation_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String connectivity__build_methylation_correlation_covariates_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String connectivity__finalize_methylation_connectivity_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String mergecovariates__merge_covariatesr_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
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
        Boolean AnnotateSites = true
        String ValueColumn = "mod_score"
        Float ValueMultiplier = 0.01
        Boolean ComputeCoverageMethylationCorrelation = true
        Int MergeMemoryGB = 128
        Int MergeDiskGB = 500
        Int AggregateMemoryGB = 64
        Int AggregateDiskGB = 1000
        Int AnnotationMemoryGB = 256
        Int AnnotationDiskGB = 200
        Int CorrelationWindowBP = 1000
        Float CorrelationMinAbsCorrelation = 0.95
        Int CorrelationMemoryGB = 64
        Int CorrelationDiskGB = 250
        Float ConnectivityZThreshold = -3.0
        Int NumThreads = 1

    }

    call CohortAggregation.AggregateMethylationData as AggregateCohort {
        input:
      annotation__annotate_methylation_sites_image = annotation__annotate_methylation_sites_image,
      cohort_aggregation__aggregate_methylation_chromosomes_image = cohort_aggregation__aggregate_methylation_chromosomes_image,
      cohort_aggregation__build_methylation_cohort_samples_image = cohort_aggregation__build_methylation_cohort_samples_image,
      cohort_aggregation__merge_methylation_chromosome_image = cohort_aggregation__merge_methylation_chromosome_image,
      cohort_aggregation__prepare_methylation_cohort_manifest_image = cohort_aggregation__prepare_methylation_cohort_manifest_image,

            CohortManifest = CohortManifest,
            OutputPrefix = OutputPrefix,
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            AutosomePrefix = AutosomePrefix,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            ComputeCoverageMethylationCorrelation = ComputeCoverageMethylationCorrelation,
            AnnotationGTF = AnnotationGTF,
            CCREAnnotations = CCREAnnotations,
            CpGIslandAnnotations = CpGIslandAnnotations,
            PromoterWindow = PromoterWindow,
            AnnotateSites = AnnotateSites,
            MergeMemoryGB = MergeMemoryGB,
            MergeDiskGB = MergeDiskGB,
            AggregateMemoryGB = AggregateMemoryGB,
            AggregateDiskGB = AggregateDiskGB,
            AnnotationMemoryGB = AnnotationMemoryGB,
            AnnotationDiskGB = AnnotationDiskGB,
            NumThreads = NumThreads
  }

    call Connectivity.RefineMethylationConnectivity as RefineConnectivity {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,
      connectivity__analyze_methylation_cpg_correlation_image = connectivity__analyze_methylation_cpg_correlation_image,
      connectivity__build_methylation_correlation_covariates_image = connectivity__build_methylation_correlation_covariates_image,
      connectivity__finalize_methylation_connectivity_image = connectivity__finalize_methylation_connectivity_image,

            IntMethylationBedsByChromosome = AggregateCohort.IntMethylationBedsByChromosome,
            ChromosomeOutputSuffixes = AggregateCohort.ChromosomeOutputSuffixes,
            PreConnectivityFilteredCalls = AggregateCohort.PreConnectivityFilteredCalls,
            PreConnectivityRawMethylationBed = AggregateCohort.PreConnectivityRawMethylationBed,
            PreConnectivityIntMethylationBed = AggregateCohort.PreConnectivityIntMethylationBed,
            PreConnectivitySampleQC = AggregateCohort.PreConnectivitySampleQC,
            AdditionalCovariates = AdditionalCovariates,
            OutputPrefix = OutputPrefix,
            PcMemoryGB = MergeMemoryGB,
            PcDiskGB = MergeDiskGB,
            NumThreads = NumThreads,
            CorrelationWindowBP = CorrelationWindowBP,
            CorrelationMinAbsCorrelation = CorrelationMinAbsCorrelation,
            CorrelationMemoryGB = CorrelationMemoryGB,
            CorrelationDiskGB = CorrelationDiskGB,
            ConnectivityZThreshold = ConnectivityZThreshold,
            ConnectivityMemoryGB = AggregateMemoryGB,
            ConnectivityDiskGB = AggregateDiskGB
  }

    call QtlCovariates.PrepareMethylationQtlCovariates as PrepareQtlCovariates {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,
      mergecovariates__merge_covariatesr_image = mergecovariates__merge_covariatesr_image,

            IntMethylationBed = RefineConnectivity.IntMethylationBed,
            AdditionalCovariates = AdditionalCovariates,
            OutputPrefix = OutputPrefix,
            PcMemoryGB = MergeMemoryGB,
            PcDiskGB = MergeDiskGB,
            NumThreads = NumThreads
  }

    output {
        File FilteredCalls = RefineConnectivity.FilteredCalls
        File SiteQC = AggregateCohort.SiteQC
        File SiteMetadata = AggregateCohort.SiteMetadata
        File PassingSiteMetadata = AggregateCohort.PassingSiteMetadata
        File SampleQC = RefineConnectivity.SampleQC
        File FilterSummary = AggregateCohort.FilterSummary
        File FilterCountsPlot = AggregateCohort.FilterCountsPlot
        File FilterUpsetPlot = AggregateCohort.FilterUpsetPlot
        File RawMethylationBed = RefineConnectivity.RawMethylationBed
        File IntMethylationBed = RefineConnectivity.IntMethylationBed
        File ConnectivityOutliers = RefineConnectivity.ConnectivityOutliers
        File ConnectivitySummary = RefineConnectivity.ConnectivitySummary
        File ConnectivityRepresentativeCpGs = RefineConnectivity.ConnectivityRepresentativeCpGs
        Array[File] CorrelationClustersByChromosome = RefineConnectivity.CorrelationClustersByChromosome
        Array[File] CorrelationSummariesByChromosome = RefineConnectivity.CorrelationSummariesByChromosome
        File? PassingSiteAnnotations = AggregateCohort.PassingSiteAnnotations
        File IntPhenotypePCsOut = PrepareQtlCovariates.IntPhenotypePCsOut
        File IntPhenotypePCsAllOut = PrepareQtlCovariates.IntPhenotypePCsAllOut
        File? IntQtlCovariates = PrepareQtlCovariates.IntQtlCovariates
        File CohortSamples = AggregateCohort.CohortSamples
        Int TotalSamples = AggregateCohort.TotalSamples
    }
}
