# Variable, global to REMP's namespace.
# This function is not exported to user space and does not need to be documented.
REMPOPTIONS <- options_manager(
  .default.AluFamily = c(
    "AluJb", "AluJo", "AluJr", "AluJr4", "AluSc", "AluSc5",
    "AluSc8", "AluSg", "AluSg4", "AluSg7", "AluSp", "AluSq", "AluSq10",
    "AluSq2", "AluSq4", "AluSx", "AluSx1", "AluSx3", "AluSx4", "AluSz",
    "AluSz6", "AluY", "AluYa5", "AluYa8", "AluYb8", "AluYb9", "AluYc",
    "AluYc3", "AluYc5", "AluYd8", "AluYf4", "AluYf5", "AluYg6", "AluYh9",
    "AluYk11", "AluYk12", "AluYk4", "FAM", "FLAM_A", "FLAM_C", "FRAM"
  ),

  .default.L1Family = c(
    "L1MC5a", "L1MB3", "L1MB5", "L1PA6", "L1P1", "L1MA8",
    "L1M5", "L1MA9", "L1PA14", "L1ME4c", "L1PA4", "L1PA7", "L1PA16", "L1PA2",
    "L1M2", "L1PREC2", "L1PB", "L1M4", "HAL1M8", "L1ME4a", "L1ME3E", "L1MA5",
    "L1MB8", "L1ME4b", "L1MC3", "L1MC", "L1M4c", "L1M2c", "L1M7", "L1MEc",
    "L1MA4", "L1MA7", "L1PA5", "L1MEd", "L1PA10", "L1PA8", "L1PA8A", "L1P4",
    "L1PA15-16", "L1ME3D", "L1MA2", "L1ME3B", "L1ME1", "L1MC4a", "L1MD3",
    "L1ME3", "L1MEg", "L1MC4", "L1M3", "L1M8", "HAL1", "L1M4b", "L1MC1",
    "L1PA3", "L1MB4", "L1MD", "L1M4a2", "L1MEf", "L1ME2", "L1PB1", "L1M3e",
    "L1P5", "L1MC5", "L1M4a1", "L1MCa", "L1ME3G", "L1MEb", "L1MB2", "L1ME3F",
    "L1MB7", "L1MCc", "L1ME3A", "L1PA12", "L1MDa", "L1ME3Cz", "L1ME2z",
    "L1PA15", "L1MA10", "L1MC2", "L1ME3C", "L1MA5A", "L1MEh", "L1PA11",
    "L1MB1", "L1MD2", "L1PA13", "L1ME5", "L1P2", "HAL1ME", "L1MA1", "L1PA17",
    "L1PB4", "L1MA6", "L1PB2", "L1P3", "L1MD1", "HAL1b", "L1MEi", "L1MEj",
    "L1PB3", "L1M6", "L1MA3", "L1MA4A", "L1M", "L1MCb", "L1P4a", "L1M1",
    "L1M3f", "L1HS", "L1M2a", "L1P4e", "L1PBa", "L1M2b", "L1MDb", "L1M3c",
    "L1M6B", "L1MEg2", "L1MEg1", "L1M3de", "L1M3d", "L1PBb", "L1M3a", "L1P3b",
    "L1P4d", "L1MEa", "L1M2a1", "L1M3b", "L1P4b", "L1PBa1", "L1P", "L1P4c"
  ),

  .default.GM12878.450k.URL = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/wgEncodeHaibMethyl450Gm12878SitesRep1.bed.gz",

  .default.AH.repeatmasker.hg19 = "AH5122",

  .default.AH.refgene.hg19 = "AH5040",

  .default.TSS.upstream = 2000,

  .default.TSS.downstream = 0,

  .default.max.flankWindow = 1200,

  .default.27k.total.probes = 27578,

  .default.450k.total.probes = 485577,

  .default.epic.total.probes = 866895,

  .default.450k.annotation = "ilmn12.hg19",

  .default.epic.annotation = "ilm10b2.hg19",

  .default.genomicRegionColNames = c(
    "InNM", "InNR", "InTSS", "In5UTR",
    "InCDS", "InExon", "In3UTR"
  ),

  .default.predictors = c(
    "RE.score", "RE.Length", "RE.CpG.density",
    "RE.InTSS", "RE.In5UTR", "RE.InCDS",
    "N.nbr", "distance.mean", "distance.std", "distance.min", "distance.min2",
    "Methy.min", "Methy.min2", "Methy.mean", "Methy.std"
  ),

  .default.svmLinear.tune = list(C = 2^seq(-15, 3, 2)),

  .default.svmRadial.tune = list(
    C = 2^seq(-7, 7, 2),
    sigma = 2^seq(-9, 1, 2)
  ),

  .default.xgbTree.tune = list(
    nrounds = 500,
    eta = 0.01,
    max_depth = 5,
    gamma = seq(0, 1, 0.2),
    colsample_bytree = 0.8,
    min_child_weight = 4,
    subsample = 0.8
  )
)


#' @title Set or get options for REMP package
#'
#' @description Tools to manage global setting options for REMP package.
#'
#' @param ... Option names to retrieve option values or \code{[key]=[value]} pairs to set options.
#'
#' @section Supported options:
#' The following options are supported
#' \describe{
#'  \item{\code{.default.AluFamily}}{A list of Alu subfamily to be included in the prediction.}
#'  \item{\code{.default.L1Family}}{A list of L1 subfamily to be included in the prediction.}
#'  \item{\code{.default.GM12878.450k.URL}}{ URL to download GM12878 450k methylation profiling data.}
#'  \item{\code{.default.AH.repeatmasker.hg19}}{\code{AnnotationHub} data ID linked to RepeatMasker
#'  database (build hg19)}
#'  \item{\code{.default.AH.refgene.hg19}}{\code{AnnotationHub} data ID linked to refSeq gene database
#'  (build hg19)}
#'  \item{\code{.default.TSS.upstream}}{Define the upstream range of transcription start site region.}
#'  \item{\code{.default.TSS.downstream}}{Define the downstream range of transcription start site region.}
#'  \item{\code{.default.max.flankWindow}}{Define the max size of the flanking window surrounding the
#'  predicted RE-CpG.}
#'  \item{\code{.default.450k.total.probes}}{Total number of probes designed in Illumina 450k array.}
#'  \item{\code{.default.epic.total.probes}}{Total number of probes designed in Illumina EPIC array.}
#'  \item{\code{.default.450k.annotation}}{A character string associated with the Illumina 450k
#'  array annotation dataset.}
#'  \item{\code{.default.epic.annotation}}{A character string associated with the Illumina EPIC
#'  array annotation dataset.}
#'  \item{\code{.default.genomicRegionColNames}}{Define the names of the genomic regions for prediction.}
#'  \item{\code{.default.predictors}}{Define the names of predictors for RE methylation prediction.}
#'  \item{\code{.default.svmLinear.tune}}{Define the default \code{C} (Cost) parameter for Support
#'  Vector Machine (SVM) using linear kernel.}
#'  \item{\code{.default.svmRadial.tune}}{Define the default parameters (\code{C} and \code{sigma}) for SVM
#'  using Radial basis function kernel.}
#'  \item{\code{.default.xgbTree.tune}}{Define the default parameters (\code{nrounds}, \code{eta}, \code{max_depth},
#'  \code{gamma}, \code{colsample_bytree}, \code{min_child_weight}, and \code{subsample}) for Extreme Gradient
#'  Boosting.}
#' }
#'
#' @return \code{NULL}
#'
#' @examples
#' # Display all default settings
#' remp_options()
#' 
#' # Display a specified setting
#' remp_options(".default.max.flankWindow")
#' 
#' # Change default maximum flanking window size to 2000
#' remp_options(.default.max.flankWindow = 2000)
#' 
#' # Reset all options
#' remp_reset()
#' @rdname options
#' @export
remp_options <- function(...) {
  # protect against the use of reserved words.
  stop_if_reserved(...)
  REMPOPTIONS(...)
}

#' @rdname options
#' @export
remp_reset <- function() {
  reset(REMPOPTIONS)
}
