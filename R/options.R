# Variable, global to REMP's namespace.
# This function is not exported to user space and does not need to be documented.
REMPOPTIONS <- settings::options_manager(
  .default.AluFamily.grep = "^Alu|^FAM|^FLAM|^FRAM",
  
  .default.L1Family.grep = "^L1|^HAL1",

  .default.LTRFamily.grep = "ERV|LTR",
  
  .default.chr = paste0("chr", c(seq_len(22), "X", "Y")),
  
  .default.GM12878.450k.URL = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/wgEncodeHaibMethyl450Gm12878SitesRep1.bed.gz",

  .default.RMSK.hg19.URL = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk.txt.gz",
  
  .default.RMSK.hg38.URL = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz",
  
  .default.refGene.hg19.URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
  
  .default.refGene.hg38.URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
  
  .default.AH.repeatmasker.hg19 = "AH5122",

  .default.AH.refgene.hg19 = "AH5040",

  .default.AH.hg38ToHg19.over.chain = "AH14108",
  
  .default.AH.hg19ToHg38.over.chain = "AH14150",
  
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
    "Methy.min", "Methy.min2", "Methy.mean.mov1", "Methy.mean.mov2",
    "Methy.mean.mov3", "Methy.mean.mov4", "Methy.std"
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
#'  \item{\code{.default.AluFamily.grep}}{Regular expression for 'grep' to extract Alu subfamily to be included in the prediction.}
#'  \item{\code{.default.L1Family.grep}}{Regular expression for 'grep' to extract L1 subfamily to be included in the prediction.}
#'  \item{\code{.default.LTRFamily.grep}}{Regular expression for 'grep' to extract Long Terminal Repeat (LTR), including Endogenous 
#'  Retrovirus (ERV) subfamily to be included in the prediction (this includes Human Endogenous Retrovirus, HERV).}
#'  \item{\code{.default.chr}}{List of human chromosome.}
#'  \item{\code{.default.GM12878.450k.URL}}{URL to download GM12878 450k methylation profiling data.}
#'  \item{\code{.default.RMSK.hg19.URL}}{URL to download RepeatMasker database in hg19 genome.}
#'  \item{\code{.default.RMSK.hg38.URL}}{URL to download RepeatMasker database in hg38 genome.}
#'  \item{\code{.default.refGene.hg19.URL}}{URL to download refSeq gene database in hg19 genome.}
#'  \item{\code{.default.refGene.hg38.URL}}{URL to download refSeq gene database in hg38 genome.}
#'  \item{\code{.default.AH.repeatmasker.hg19}}{\code{AnnotationHub} data ID linked to RepeatMasker
#'  database (build hg19).}
#'  \item{\code{.default.AH.refgene.hg19}}{\code{AnnotationHub} data ID linked to refSeq gene database
#'  (build hg19)}
#'  \item{\code{.default.AH.hg38ToHg19.over.chain}}{\code{AnnotationHub} hg38 to hg19 liftover chain data ID.}
#'  \item{\code{.default.AH.hg19ToHg38.over.chain}}{\code{AnnotationHub} hg19 to hg38 liftover chain data ID.}
#'  \item{\code{.default.TSS.upstream}}{Define the upstream range of transcription start site region.}
#'  \item{\code{.default.TSS.downstream}}{Define the downstream range of transcription start site region.}
#'  \item{\code{.default.max.flankWindow}}{Define the max size of the flanking window surrounding the
#'  predicted RE-CpG.}
#'  \item{\code{.default.27k.total.probes}}{Total number of probes designed in Illumina 27k array.}
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
  settings::stop_if_reserved(...)
  REMPOPTIONS(...)
}

#' @rdname options
#' @export
remp_reset <- function() {
  settings::reset(REMPOPTIONS)
}
