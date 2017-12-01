#' @title RE Annotation Database Initialization
#'
#' @description
#' \code{initREMP} is used to initialize annotation database for RE methylation prediction. 
#' Two major RE types in human, Alu element (Alu) and LINE-1 (L1) are available.
#'
#' @param arrayType Illumina methylation array type. Currently \code{"450k"} and \code{"EPIC"} 
#' are supported. Default = \code{"450k"}.
#' @param REtype Type of RE. Currently \code{"Alu"} and \code{"L1"} are supported.
#' @param RE A \code{\link{GRanges}} object containing user-specified RE genomic location information. 
#' If \code{NULL}, the function will retrive RepeatMasker RE database from \code{\link{AnnotationHub}} 
#' (build hg19).
#' @param ncore Number of cores to run parallel computation. By default max number of cores 
#' available in the machine will be utilized. If \code{ncore = 1}, no parallel computation is allowed.
#' @param BPPARAM An optional \code{\link{BiocParallelParam}} instance determining the parallel back-end to 
#' be used during evaluation. If not specified, default back-end in the machine will be used.
#' @param export Logical. Should the returned \code{\link{REMParcel}} object be saved to local machine?
#' See Details.
#' @param work.dir Path to the directory where the generated data will be saved. Valid when 
#' \code{export = TRUE}. If not specified and \code{export = TRUE}, temporary directory \code{tempdir()} 
#' will be used. 
#' @param verbose Logical parameter. Should the function be verbose?
#' 
#' @details
#' Currently, we support two major types of RE in human, Alu and L1. The main purpose of 
#' \code{initREMP} is to generate and annotate CpG/RE data using the refSeq Gene 
#' annotation database (provided by \code{\link{AnnotationHub}}). These annotation data are crucial to 
#' RE methylation prediction in \code{\link{remp}}. Once generated, the data can be reused in the future 
#' (data can be very large). Therefore, we recommend user to save the output from 
#' \code{initREMP} to the local machine, so that user only need to run this function once 
#' as long as there is no change to the RE database. To minimize the size of resulting data file, the generated 
#' annotation data are only for REs that contain RE-CpGs with neighboring profiled CpGs. By default, the 
#' neighboring CpGs are confined within 1200 bp flanking window. This window size can be modified using 
#' \code{\link{remp_options}}.
#'
#' @return An \code{\link{REMParcel}} object containing data needed for RE methylation prediction.
#'
#' @seealso See \code{\link{remp}} for RE methylation prediction.
#'
#' @examples
#' data(Alu.demo)
#' remparcel <- initREMP(arrayType = '450k', REtype = 'Alu', RE = Alu.demo, ncore = 1)
#' remparcel
#' 
#' # Save the data for later use
#' saveParcel(remparcel)
#' 
#' @export
initREMP <- function(arrayType = c("450k", "EPIC"), REtype = c("Alu", "L1"), RE = NULL, 
                     ncore = NULL, BPPARAM = NULL, 
                     export = FALSE, work.dir = tempdir(), 
                     verbose = FALSE) {
  ## Initiate running time
  t <- Sys.time()
  
  arrayType <- match.arg(arrayType)
  REtype <- match.arg(REtype)
  
  message("Start ", REtype, " annotation data initialization ...", .timeTrace(t)) 
  message("Illumina platform: ", arrayType)
  
  if (is.null(ncore)) 
    ncore <- parallel::detectCores()
  
  ## Setup backend for paralell computing
  be <- getBackend(ncore, BPPARAM, verbose)
  
  if (arrayType == "450k") {
    if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
      ILMN.GR <- minfi::getLocations(
        IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
      }
  } else if (arrayType == "EPIC") {
    if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", quietly = TRUE)) {
      ILMN.GR <- minfi::getLocations(
        IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      }
  } else stop("Wrong Illumina platform type. Can be either '450k' or 'EPIC'.")
  ILMN.GR <- ILMN.GR[substring(names(ILMN.GR), 1, 2) != "ch"]
  ILMN.GR$Index <- names(ILMN.GR)
  
  if (verbose) 
    message("Valid RE-CpG flanking window size: ", 
            remp_options(".default.max.flankWindow"), " bp")
  
  ################################################## 
  ###------ Part I - Resources preparation ------###
  ################################################## 
  
  ## All database used are from AnnotationHub(), if user does not specify
  ## the data.
  
  ## Test the file permission
  permission <- file.access(AnnotationHub::getAnnotationHubOption("CACHE"),2)
  
  if(permission != 0)
  {
    if (verbose) message(AnnotationHub::getAnnotationHubOption("CACHE"), 
                         " is not writable, using temporal directory ", 
                         .forwardSlashPath(tempdir()), " instead.")
    AnnotationHub::setAnnotationHubOption("CACHE", file.path(tempdir(),".AnnotationHub"))
  }
  
  ah <- suppressMessages(AnnotationHub::AnnotationHub())
  
  if (is.null(RE)) {
    ### Get RE annotation database (RepeatMasker)
    RE.hg19 <- fetchRMSK(ah, REtype, verbose)
  } else {
    .isGROrStop(RE)
    RE.hg19 <- RE
  }
  
  ### Get refSeq gene database
  refgene.hg19 <- fetchRefSeqGene(ah, mainOnly = FALSE, verbose)
  
  ### Get RE-CpG location database Narrow down RE.hg19 to RE sequence that
  ### overlaps with CpG sites flanking region For demo data, this will not
  ### change anything.
  ILMN.GR.flank <- .twoWayFlank(ILMN.GR, remp_options(".default.max.flankWindow"))
  RE.hg19 <- subsetByOverlaps(RE.hg19, ILMN.GR.flank, ignore.strand = TRUE)
  
  ## Narrow down RE.hg19 to RE sequence that overlaps with gene sequence
  ## RE.hg19 <- subsetByOverlaps(RE.hg19, refgene.hg19$main, ignore.strand
  ## = TRUE)
  
  ## Locate RE-CpG
  RE.CpG <- findRECpG(RE.hg19, REtype, be, verbose)
  
  ## Make RE.hg19 and RE.CpG dataset consistent with RE
  RE.hg19 <- RE.hg19[runValue(match(RE.CpG$Index, RE.hg19$Index))]
  
  ######################################## 
  ###------ Part II - Annotation ------###
  ######################################## 
  
  ### Annotate ILMN CpG
  ## Narrow down Illumina CpG Probe using RE.CpG flanking region
  RE.CpG.flanking <- .twoWayFlank(RE.CpG, remp_options(".default.max.flankWindow"))
  ILMN.GR <- subsetByOverlaps(ILMN.GR, RE.CpG.flanking, ignore.strand = TRUE)

  ###----------------------------------------------------------------------------------------------------
  ### Annotate RE
  RE.refGene <- GRannot(RE.hg19, refgene.hg19, symbol = FALSE, verbose = verbose)
  
  ###----------------------------------------------------------------------------------------------------
  ### RE-CpG covered by ILMN
  RECpG_Platform.hits <- findOverlaps(ILMN.GR, RE.CpG, ignore.strand = TRUE)
  mcols(ILMN.GR)$RE.Index <- Rle(NA)
  ILMN.GR$RE.Index[queryHits(RECpG_Platform.hits)] <- RE.CpG[subjectHits(RECpG_Platform.hits)]$Index
  
  remparcel <- REMParcel(REtype = REtype, platform = arrayType,
                         RefGene = refgene.hg19$main, 
                         RE = RE.refGene, RECpG = RE.CpG,
                         ILMN = ILMN.GR)
  
  message("Done.", .timeTrace(t))
  
  if(export)
    saveParcel(remparcel, work.dir, verbose)
  
  return(remparcel)
}  # End of intREMP
