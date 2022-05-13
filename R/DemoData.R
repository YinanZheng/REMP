#' Subset of Alu genomic location dataset (hg19)
#'
#' A \code{\link{GRanges}} dataset containing 500 Alu sequences that have CpGs profiled
#' by both Illumina 450k and EPIC array. The variables are as follows:
#'
#' \itemize{
#'   \item seqnames: chromosome number
#'   \item ranges: hg19 genomic position
#'   \item strand: DNA strand
#'   \item swScore: Smith Waterman (SW) alignment score
#'   \item repName: Alu name
#'   \item repClass: Alu class
#'   \item repFamily: Alu family
#'   \item Index: internal index (meaningless for external use; 
#'   not communicable between genome builds)
#' }
#'
#' \code{Alu.hg19.demo} has the same format as the data object 
#' returned by \code{\link{fetchRMSK}}.
#'
#' @format A \code{\link{GRanges}} object.
#' 
#' @return A GRanges object with 500 ranges and 3 metadata columns.
#' 
#' @seealso See \code{\link{fetchRMSK}} to obtain the complete Alu/L1 dataset.

#' @source RepeatMasker database provided by package \code{AnnotationHub}
"Alu.hg19.demo"

#' Subset of Alu genomic location dataset (hg38)
#'
#' A \code{\link{GRanges}} dataset containing 500 Alu sequences that have CpGs profiled
#' by both Illumina 450k and EPIC array. The variables are as follows:
#'
#' \itemize{
#'   \item seqnames: chromosome number
#'   \item ranges: hg38 genomic position
#'   \item strand: DNA strand
#'   \item swScore: Smith Waterman (SW) alignment score
#'   \item repName: Alu name
#'   \item repClass: Alu class
#'   \item repFamily: Alu family
#'   \item Index: internal index (meaningless for external use; 
#'   not communicable between genome builds)
#' }
#'
#' \code{Alu.hg38.demo} has the same format as the data object 
#' returned by \code{\link{fetchRMSK}}.
#'
#' @format A \code{\link{GRanges}} object.
#' 
#' @return A GRanges object with 500 ranges and 3 metadata columns.
#' 
#' @seealso See \code{\link{fetchRMSK}} to obtain the complete Alu/L1 dataset.

#' @source RepeatMasker database (hg38) is provided by UCSC database downloaded 
#' from http://hgdownload.cse.ucsc.edu/goldenpath.
"Alu.hg38.demo"
