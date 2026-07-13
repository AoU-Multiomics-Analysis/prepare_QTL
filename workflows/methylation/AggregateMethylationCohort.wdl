version 1.0
import "cohort_aggregation.wdl" as CohortAggregation
import "connectivity.wdl" as Connectivity
import "qtl_covariates.wdl" as QtlCovariates

# Public cohort-level entry point for per-sample methylation outputs. The input
# and output contract stays stable while implementation is split into logical
# stage workflows in this directory.

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

    call CohortAggregation.AggregateMethylationData as AggregateCohort {
        input:
            CohortManifest = CohortManifest,
            OutputPrefix = OutputPrefix,
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            AutosomePrefix = AutosomePrefix,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            AnnotationGTF = AnnotationGTF,
            CCREAnnotations = CCREAnnotations,
            CpGIslandAnnotations = CpGIslandAnnotations,
            PromoterWindow = PromoterWindow,
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
            MaxConnectivityFeatures = MaxConnectivityFeatures,
            ConnectivityLandmarks = ConnectivityLandmarks,
            ConnectivityZThreshold = ConnectivityZThreshold,
            ConnectivityMemoryGB = AggregateMemoryGB,
            ConnectivityDiskGB = AggregateDiskGB
    }

    call QtlCovariates.PrepareMethylationQtlCovariates as PrepareQtlCovariates {
        input:
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
        File PassingSiteAnnotations = AggregateCohort.PassingSiteAnnotations
        File IntPhenotypePCsOut = PrepareQtlCovariates.IntPhenotypePCsOut
        File? IntQtlCovariates = PrepareQtlCovariates.IntQtlCovariates
        File CohortSamples = AggregateCohort.CohortSamples
        Int TotalSamples = AggregateCohort.TotalSamples
    }
}
