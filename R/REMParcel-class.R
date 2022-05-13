#' @title REMParcel instances
#'
#' @description \code{REMParcel} is a container class to organize required datasets for
#' RE methylation prediction generated from \code{\link{initREMP}} and used in \code{\link{remp}}.
#'
#' @name REMParcel-class
#'
#' @rdname REMParcel-class
#'
#' @param object A \code{REMParcel} object.
#' @param REtype Type of RE (\code{"Alu"}, \code{"L1"}, or \code{"ERV"}).
#' @param genome Specify the build of human genome. Can be either \code{"hg19"} or \code{"hg38"}.
#' @param platform Illumina methylation profiling platform (\code{"450k"} or \code{"EPIC"}).
#' @param RefGene refSeq gene annotation data, which can be obtained by \code{\link{fetchRefSeqGene}}.
#' @param RE Annotated RE genomic range data, which can be obtained by \code{\link{fetchRMSK}} and annotated
#' by \code{\link{GRannot}}.
#' @param RECpG Genomic range data of annotated CpG site identified in RE DNA sequence, which can
#' be obtained by \code{\link{findRECpG}} and annotated by \code{\link{GRannot}}.
#' @param ILMN Illumina CpG probe genomic range data.
#' @param REonly For \code{getILMN}: see Accessors.
#' @param work.dir For \code{saveParcel}: path to the directory where the generated data will be saved.
#' If not specified, temporary directory \code{tempdir()} will be used.
#' @param verbose For \code{saveParcel}: logical parameter. Should the function be verbose?
#' @param ... For \code{saveParcel}: other parameters to be passed to the \code{saveRDS} method. See
#' \code{\link{saveRDS}}.
#'
#' @return An object of class \code{REMParcel} for the constructor.
#'
#' @section Accessors:
#' \describe{
#'     \item{\code{getParcelInfo(object)}}{Return data type, RE type, and flanking window size information
#'     of the parcel.}
#'     \item{\code{getRefGene(object)}}{Return RefSeq gene annotation data.}
#'     \item{\code{getRE(object)}}{Return RE genomic location data for prediction
#'     (annotated by refSeq gene database).}
#'     \item{\code{getRECpG(object)}}{Return RE-CpG genomic location
#'     data for prediction.}
#'     \item{\code{getILMN(object, REonly = FALSE)}}{Return Illumina CpG probe genomic
#'     location data for prediction (annotated by refSeq gene database). If
#'     \code{REonly = TRUE}, only probes within RE region are returned.}
#'     }
#'
#' @section Utilities:
#' \describe{
#'     \item{\code{saveParcel(object, work.dir = tempdir(), verbose = FALSE, ...)}}{Save
#'     the object to local machine. }
#'     }
#'
#' @examples
#' showClass("REMParcel")
#' @exportClass REMParcel
REMParcel <- setClass("REMParcel",
  slots = c(
    REMParcelInfo = "CharacterList",
    RefGene = "GRanges",
    RE = "GRanges",
    RECpG = "GRanges",
    ILMN = "GRanges"
  ),
  validity = function(object) {
    if (!identical(
      runValue(object@RE$Index),
      runValue(object@RECpG$Index)
    )) {
      stop("RE data does not match RE-CpG data.")
    }
  }
)

## Constructor function
REMParcel <- function(REtype = "Unknown", 
                      genome = "Unknown", 
                      platform = "Unknown",
                      RefGene = GRanges(),
                      RE = GRanges(), RECpG = GRanges(), ILMN = GRanges()) {
  REMParcelInfo <- CharacterList(
    REtype = REtype,
    genome = genome,
    platform = platform,
    max.win = as.character(remp_options(".default.max.flankWindow"))
  )

  new("REMParcel",
    REMParcelInfo = REMParcelInfo,
    RefGene = RefGene,
    RE = RE,
    RECpG = RECpG,
    ILMN = ILMN
  )
}
