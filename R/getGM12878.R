#' @title Get methylation data of HapMap LCL sample GM12878 profiled by 
#' Illumina 450k array or EPIC array
#'
#' @description
#' \code{getGM12878} is used to obtain public available methylation profiling 
#' data of HapMap LCL sample GM12878.
#'
#' @param arrayType Illumina methylation array type. Currently \code{"450k"} 
#' and \code{"EPIC"} are supported. Default = \code{"450k"}.
#' @param mapGenome logical parameter. If \code{TRUE}, function will return 
#' a \code{\link{GenomicRatioSet}} object instead of an \code{link{RatioSet}} 
#' object.

#' @details
#' Illumina 450k data were sourced and curated from ENCODE 
#' \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/wgEncodeHaibMethyl450Gm12878SitesRep1.bed.gz}. 
#' Illumina EPIC data were obtained from data package \code{minfiDataEPIC}.
#' 
#' @return An \code{\link{RatioSet}} or \code{\link{GenomicRatioSet}} containing 
#' beta value and M value of the methylation data.
#'
#' @examples
#' # Get GM12878 methylation data (EPIC array)
#' GM12878_EPIC <- getGM12878('EPIC') 
#' GM12878_EPIC
#' 
#' @export
getGM12878 <- function(arrayType = c("450k", "EPIC"), mapGenome = FALSE) {
  arrayType <- match.arg(arrayType)
  
  if (toupper(arrayType) == "450K") {
    ### GM12878-450k
    tempfilenames <- list.files(tempdir())
    which.GM12878.450k.file <- grep("remp.GM12878.450k.file", tempfilenames)
    if (length(which.GM12878.450k.file) == 0) {
      temp <- tempfile(pattern = "remp.GM12878.450k.file")  # Make an identifiable ext to avoid repeated download
      download.file(remp_options(".default.GM12878.450k.URL"), temp, quiet = TRUE)
    } else {
      temp <- file.path(tempdir(), tempfilenames[which.GM12878.450k.file])
    }
    GM12878_450k <- read.table(gzfile(temp), header = FALSE)
    GM12878_450k.d <- matrix(GM12878_450k[, "V5"]/1000, ncol = 1)
    colnames(GM12878_450k.d) <- "GM12878"
    rownames(GM12878_450k.d) <- GM12878_450k[, "V4"]
    GM12878.RatioSet <- minfi::RatioSet(Beta = GM12878_450k.d, M = .toM(GM12878_450k.d), 
                                        annotation = c(array = "IlluminaHumanMethylation450k", 
                                                       annotation = remp_options(".default.450k.annotation")))
    if (mapGenome) 
      GM12878.RatioSet <- minfi::mapToGenome(GM12878.RatioSet)
  }
  
  if (toupper(arrayType) == "EPIC") {
    ### GM12878-850k
    requireNamespace("minfiDataEPIC", quietly = TRUE)
    GM12878_EPIC.beta <- minfi::getBeta(minfiDataEPIC::MsetEPIC)

    GM12878_EPIC.beta <- rowMeans(GM12878_EPIC.beta, na.rm = TRUE)
    GM12878_EPIC.M <- .toM(GM12878_EPIC.beta)
    
    GM12878_EPIC.beta.d <- matrix(GM12878_EPIC.beta, ncol = 1)
    GM12878_EPIC.M.d <- matrix(GM12878_EPIC.M, ncol = 1)
    
    colnames(GM12878_EPIC.beta.d) <- "GM12878"
    colnames(GM12878_EPIC.M.d) <- "GM12878"
    
    rownames(GM12878_EPIC.beta.d) <- names(GM12878_EPIC.beta)
    rownames(GM12878_EPIC.M.d) <- names(GM12878_EPIC.M)
    
    GM12878.RatioSet <- minfi::RatioSet(Beta = GM12878_EPIC.beta.d, 
                                        M = GM12878_EPIC.M.d, 
                                        annotation = c(array = "IlluminaHumanMethylationEPIC", 
                                                       annotation = remp_options(".default.epic.annotation")))
    if (mapGenome) 
      GM12878.RatioSet <- minfi::mapToGenome(GM12878.RatioSet)
  }
  
  return(GM12878.RatioSet)
}  ## End of getGM12878