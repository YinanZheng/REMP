#' Subset of Alu sequence dataset
#'
#' A \code{GRanges} dataset containing 500 Alu sequences that have CpGs profiled 
#' by both Illumina 450k and EPIC array. The variables are as follows:
#'
#' \itemize{
#'   \item seqnames chromosome number
#'   \item ranges hg19 genomic position
#'   \item strand DNA strand
#'   \item name Alu subfamily
#'   \item score Smith Waterman (SW) alignment score
#'   \item Index internal index
#' }
#'
#' Scripts for generating this object are contained in the \code{scripts} directory.
#'
#' @format A \code{GRanges} object.
#' @return A GRanges object with 500 ranges and 3 metadata columns.
#' @source RepeatMasker database provided by package \code{AnnotationHub}.
"Alu.demo"
