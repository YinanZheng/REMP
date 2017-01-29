#' @title Repetitive Element Methylation Prediction
#'
#' @description
#' Machine learing-based tools to predict DNA methylation of locus-specific repetitive elements (RE) 
#' by learning surrounding genetic and epigenetic information. These tools provide genomewide and 
#' single-base resolution of DNA methylation prediction on RE that are difficult to measure using 
#' array-based or sequencing-based platforms, which enables epigenome-wide association study (EWAS) 
#' and differentially methylated region (DMR) analysis on RE.
#' 
#' \tabular{ll}{ Name: \tab REMP\cr Type: \tab Package\cr
#' Version: \tab 0.99.16\cr License: \tab GPL-3\cr }
#' 
#' @name REMP-package
#' 
#' @aliases REMP-package REMP
#' 
#' @docType package
#' 
#' @section Overview - standard procedure: 
#' \describe{
#'     \item{Step 1}{Start out generating required dataset for prediction using \code{\link{initREMP}}. 
#'     It is recommended to save these generated data to working directory so they can be used in the future.}
#'     \item{Step 2}{Clean Illumina methylation dataset using \code{\link{grooMethy}}. This function 
#'     can help identify and fix abnormal values and automatically impute missing values, which are 
#'     essential for downstream prediction.}
#'     \item{Step 3}{Run \code{\link{remp}} to predict genome-wide locus specific RE methylation.}
#'     \item{Step 4}{Use the built-in accessors and utilities in \code{\link{REMProduct}} object to get or 
#'     refine the prediction results.}
#'     }
#'     
#' @author 
#' Yinan Zheng \email{y-zheng@@northwestern.edu},
#' Lei Liu \email{lei.liu@@northwestern.edu},       
#' Wei Zhang \email{wei.zhang1@@northwestern.edu},
#' Warren Kibbe \email{warren.kibbe@@nih.gov},       
#' Lifang Hou \email{l-hou@@northwestern.edu}
#' 
#' Maintainer: Yinan Zheng \email{y-zheng@@northwestern.edu}
#' 
#' @references Manuscript in review (Nucleic Acids Research)
#' 
#' @keywords package
#' 
#' @import methods
#' @import S4Vectors
#' @import BiocParallel 
#' @import caret 
#' @import impute 
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' 
#' @importFrom quantregForest quantregForest
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges assays rowData colData
#' @importFrom settings stop_if_reserved reset options_manager
#' @importFrom stats setNames density predict
#' @importFrom utils download.file read.table data
#' @importFrom graphics lines plot
#' @importFrom AnnotationHub AnnotationHub getAnnotationHubOption setAnnotationHubOption
#' @importFrom iterators nextElem idiv
#' @importFrom org.Hs.eg.db org.Hs.egREFSEQ2EG org.Hs.egSYMBOL
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom GenomicRanges GRanges GRangesList makeGRangesListFromFeatureFragments findOverlaps
#' @importFrom GenomicRanges granges seqnames start end strand ranges promoters shift start<- end<-
#' @importFrom IRanges IRanges IRangesList CharacterList subsetByOverlaps
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom BiocGenerics mget
#' @importFrom BSgenome getSeq 
#' @importFrom minfi getAnnotation getBeta getM MethylSet RatioSet GenomicRatioSet IlluminaMethylationAnnotation
#' @importFrom Biostrings DNAString vmatchPattern

NULL