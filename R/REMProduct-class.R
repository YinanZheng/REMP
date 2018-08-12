#' @title REMProduct instances
#' 
#' @description Class \code{REMProduct} is to maintain RE methylation prediction results. 
#' \code{REMProduct} inherits Bioconductor's \code{RangedSummarizedExperiment} class.
#' 
#' @name REMProduct-class
#' 
#' @rdname REMProduct-class
#' 
#' @param object A \code{REMProduct} object.
#' @param object1 A \code{REMProduct} object.
#' @param object2 A \code{REMProduct} object.
#' @param REtype Type of RE (\code{"Alu"} or \code{"L1"}).
#' @param platform Illumina methylation profiling platform (\code{"450k"} or \code{"EPIC"}).
#' @param win Flanking window size of the predicting RE-CpG.
#' @param predictModel Name of the model used for prediction.
#' @param QCModel Name of the model used for prediction quality evaluation.
#' @param rempM Predicted methylation level in M value.
#' @param rempB Predicted methylation level in beta value (optional).
#' @param rempQC Prediction quality scores, which is available only when Random Forest 
#' model is used in \code{\link{remp}}.
#' @param cpgRanges Genomic ranges of the predicting RE-CpG.
#' @param sampleInfo Sample information. 
#' @param REannotation Annotation data for the predicting RE.
#' @param RECpG Annotation data for the RE-CpG profiled by Illumina platform.
#' @param regionCode Internal index code defined in \code{refGene} for gene region indicators.
#' @param refGene refSeq gene annotation data, which can be obtained by \code{\link{fetchRefSeqGene}}.
#' @param varImp Importance of the predictors.
#' @param REStats RE coverage statistics, which is internally generated in \code{\link{remp}}.
#' @param GeneStats Gene coverage statistics, which is internally generated in \code{\link{remp}}.
#' @param Seed Random seed for Random Forest model for reproducible prediction results.
#' @param type For \code{plot} and \code{decodeAnnot}: see Utilities.
#' @param ncore For \code{decodeAnnot} and \code{rempAggregate}: number of cores to run parallel computation. 
#' By default no parallel computation is allowed (\code{ncore = 1}).
#' @param BPPARAM For \code{decodeAnnot} and \code{rempAggregate}: an optional \code{\link{BiocParallelParam}} 
#' instance determining the parallel back-end to be used during evaluation. If not specified, default 
#' back-end in the machine will be used.
#' @param x For \code{plot}: an \code{REMProduct} object.
#' @param threshold For \code{rempTrim}: see Utilities.
#' @param missingRate For \code{rempTrim}: see Utilities.
#' @param NCpG For \code{rempAggregate}: see Utilities.
#' @param ... For \code{plot}: \code{\link{graphical parameters}} to be passed to the \code{plot} method.
#' 
#' @return An object of class \code{REMProduct} for the constructor.
#' 
#' @section Accessors:
#' \describe{
#'     \item{\code{rempM(object)}}{Return M value of the prediction.}
#'     \item{\code{rempB(object)}}{Return beta value of the prediction.}
#'     \item{\code{rempQC(object)}}{Return prediction quality scores.}
#'     \item{\code{rempImp(object)}}{Return relative importance of predictors.}
#'     \item{\code{rempStats(object)}}{Return RE and gene coverage statistics.}
#'     \item{\code{rempAnnot(object)}}{Return annotation data for the predicted RE.}     
#'     }
#'     
#' @section Utilities:
#' \describe{
#'     \item{\code{plot(x, type = c("individual", "overall"), ...)}}{Make a density plot of predicted methylation 
#'     (beta values) in the \code{REMProduct} object \code{x}. If \code{type = "individual"}, density curves will be 
#'     plotted for each of the samples; If \code{type = "overall"}, one density curve of the mean methylation level 
#'     across the samples will be plotted. Default \code{type = "individual"}.}
#'     \item{\code{details(object)}}{Display detailed descriptive statistics of the predicion results.}
#'     \item{\code{decodeAnnot(object, type = c("symbol", "entrez")), ncore = NULL, BPPARAM = NULL}}{Decode the 
#'     RE annotation data by Gene Symbol (when \code{type = "Symbol"}) or Entrez Gene 
#'     (when \code{type = "Entrez"}).Default \code{type = "Symbol"}. Annotation data are provided by 
#'     \code{\link{org.Hs.eg.db}}.}
#'     \item{\code{rempTrim(object, threshold = 1.7, missingRate = 0.2)}}{Any predicted CpG values with 
#'     quality score < threshold (default = 1.7, specified by \code{threshold = 1.7}) will be replaced with NA. 
#'     CpGs contain more than missingRate * 100% (default = 20%, specified by \code{missingRage = 0.2}) missing 
#'     rate across samples will be discarded. Relavant statistics will be re-evaluated.}
#'     \item{\code{rempAggregate(object, NCpG = 2, ncore = NULL, BPPARAM = NULL)}}{Aggregate the predicted RE-CpG 
#'     methylation by RE using mean. To ensure the reliability of the aggregation, by default only RE with at 
#'     least 2 predicted CpG sites (specified by \code{NCpG = 2}) will be aggregated.}
#'     \item{\code{rempCombine(object1, object2)}}{Combine two \code{REMProduct} objects by column.}
#'     }
#' 
#' @examples
#' showClass("REMProduct")
#' 
#' @exportClass REMProduct
REMProduct <- setClass("REMProduct", 
                    slots = c(REMPInfo = "CharacterList"),
                    contains = "RangedSummarizedExperiment"
                    )

## Constructor function
REMProduct <- function(REtype = "Unknown", platform = "Unknown", win = "Unknown", 
                       predictModel = "Unknown", QCModel = "Unknown",
                       rempM = NULL, rempB = NULL, rempQC = NULL, 
                       cpgRanges = GRanges(), sampleInfo = DataFrame(),
                       REannotation = GRanges(), RECpG = GRanges(), 
                       regionCode = DataFrame(), 
                       refGene = GRanges(),
                       varImp = DataFrame(),
                       REStats = DataFrame(), GeneStats = DataFrame(),
                       Seed = NULL)
{
  rempInfo <- CharacterList(REtype = REtype,
                            platform = platform,
                            win = win, 
                            predictModel = predictModel,
                            QCModel = QCModel)
  if(is.null(rempB)) rempB <- .toBeta(rempM)
  assays <- SimpleList(rempB = rempB, rempM = rempM, rempQC = rempQC)
  assays <- assays[!vapply(assays, is.null, logical(1))]
  
  stopifnot(identical(REannotation, unique(REannotation)))
  
  metadata = list(REannotation = REannotation,
                  RECpG = RECpG,
                  regionCode = regionCode,
                  refGene = refGene,
                  varImp = varImp,
                  REStats = REStats,
                  GeneStats = GeneStats,
                  Seed = Seed)
  new("REMProduct",
      REMPInfo = rempInfo,
      SummarizedExperiment(
        assays = assays,
        rowRanges = cpgRanges,
        colData = sampleInfo,
        metadata = metadata)
  )
}
