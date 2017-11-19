#' @title Repetitive Element Methylation Prediction
#'
#' @description
#' Machine learning-based tools to predict DNA methylation of locus-specific repetitive elements (RE) 
#' by learning surrounding genetic and epigenetic information. These tools provide genomewide and 
#' single-base resolution of DNA methylation prediction on RE that are difficult to measure using 
#' array-based or sequencing-based platforms, which enables epigenome-wide association study (EWAS) 
#' and differentially methylated region (DMR) analysis on RE.
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
#'     The datasets include RE information, RE-CpG (i.e. CpGs located in RE region) information, 
#'     and gene annotation, which are maintained in a \code{\link{REMParcel}} object. 
#'     It is recommended to save these generated data to the working directory so they can be used in the future.}
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
#' @references Zheng Y, Joyce BT, Liu L, Zhang Z, Kibbe WA, Zhang W, Hou L. 
#' Prediction of genome-wide DNA methylation in repetitive elements. 
#' Nucleic Acids Res. 2017;45(15):8697-711. 
#' PubMed PMID: 28911103; PMCID: PMC5587781.
#' http://dx.doi.org/10.1093/nar/gkx587.
#' 
#' @keywords package
#' 
#' @import methods
#' @import S4Vectors
#' @import BiocParallel 
#' @import caret
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges assays rowData colData
#' @importFrom minfi getAnnotation getBeta getM MethylSet RatioSet GenomicRatioSet IlluminaMethylationAnnotation
#' @importFrom impute impute.knn
#' @importFrom kernlab predict
#' @importFrom ranger ranger importance
#' @importFrom settings stop_if_reserved reset options_manager
#' @importFrom stats setNames density na.pass
#' @importFrom utils download.file read.table data
#' @importFrom graphics lines plot
#' @importFrom AnnotationHub AnnotationHub getAnnotationHubOption setAnnotationHubOption
#' @importFrom iterators nextElem idiv
#' @importFrom org.Hs.eg.db org.Hs.egREFSEQ2EG org.Hs.egSYMBOL
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 Hsapiens
#' @importFrom GenomicRanges GRanges GRangesList makeGRangesListFromFeatureFragments findOverlaps
#' @importFrom GenomicRanges granges seqnames start end strand ranges promoters shift start<- end<-
#' @importFrom IRanges IRanges IRangesList CharacterList NumericList subsetByOverlaps
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom BiocGenerics mget
#' @importFrom BSgenome getSeq 
#' @importFrom Biostrings DNAString vmatchPattern

NULL