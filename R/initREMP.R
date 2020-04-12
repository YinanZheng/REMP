#' @title RE Annotation Database Initialization
#'
#' @description
#' \code{initREMP} is used to initialize annotation database for RE methylation prediction.
#' Two major RE types in human, Alu element (Alu) and LINE-1 (L1) are available.
#'
#' @param arrayType Illumina methylation array type. Currently \code{"450k"}, \code{"EPIC"},
#' and \code{"Sequencing"} are supported. Default = \code{"450k"}.
#' @param REtype Type of RE. Currently \code{"Alu"}, \code{"L1"}, and \code{"LTR"} are supported.
#' @param annotation.source Character parameter. Specify the source of annotation databases, including
#' the RefSeq Gene annotation database and RepeatMasker annotation database. If \code{"AH"}, the database 
#' will be obtained from the AnnotationHub package. If \code{"UCSC"}, the database will be downloaded 
#' from the UCSC website http://hgdownload.cse.ucsc.edu/goldenpath. The corresponding build (\code{"hg19"} or 
#' \code{"hg38"}) can be specified in the parameter \code{genome}.
#' @param genome Character parameter. Specify the build of human genome. Can be either \code{"hg19"} or 
#' \code{"hg38"}. Note that if \code{annotation.source == "AH"}, only hg19 database is available.
#' @param RE A \code{\link{GRanges}} object containing user-specified RE genomic location information.
#' If \code{NULL}, the function will retrive RepeatMasker RE database from \code{\link{AnnotationHub}}
#' (build hg19) or download the database from UCSC website (build hg19/hg38).
#' @param Seq.GR A \code{\link{GRanges}} object containing genomic locations of the CpGs profiled by sequencing
#' platforms. This parameter should not be \code{NULL} if \code{arrayType == 'Sequencing'}. Note that the genomic
#' location can be in either hg19 or hg38 build. See details.
#' @param ncore Number of cores used for parallel computing. By default max number of cores
#' available in the machine will be utilized. If \code{ncore = 1}, no parallel computing is allowed.
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
#' Currently, we support two major types of RE in the human genome, Alu and L1. The main purpose of
#' \code{initREMP} is to generate and annotate CpG/RE data using the refSeq Gene (hg19)
#' annotation database (provided by \code{\link{AnnotationHub}}). These annotation data are crucial to
#' RE methylation prediction in \code{\link{remp}}. Once generated, the data can be reused in the future
#' (data can be very large). Therefore, we recommend the user to save the output from
#' \code{initREMP} to the local machine, so that user only need to run this function once
#' as long as there is no change to the RE database. To minimize the size of the resulting data file, the generated
#' annotation data are only for REs that contain RE-CpGs with neighboring profiled CpGs. By default, the
#' neighboring CpGs are confined within 1200 bp flanking window. This window size can be modified using
#' \code{\link{remp_options}}. Note that the refSeq Gene database from UCSC is dynamic (updated periodically) 
#' and reflecting the latest knowledge of gene, whereas the database from AnnotationHub is static and classic. 
#' Using different sources will have a slight impact on the prediction results of RE methylation and gene annotation 
#' of final results. For sequencing methylation data, please specify the genomic location of CpGs
#' in a \code{GenomicRanges} object and specify it in \code{Seq.GR}. For an example of \code{Seq.GR}, Please
#' run \code{minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)} (the row names of the CpGs in
#' \code{Seq.GR} can be \code{NULL}). The user should make sure the genome build of \code{Seq.GR} match the 
#' build specified in \code{genome} parameter (default is \code{"hg19"}).
#'
#' @return An \code{\link{REMParcel}} object containing data needed for RE methylation prediction.
#'
#' @seealso See \code{\link{remp}} for RE methylation prediction.
#'
#' @examples
#' if (!exists("remparcel")) {
#'   data(Alu.hg19.demo)
#'   remparcel <- initREMP(arrayType = "450k", 
#'                         REtype = "Alu", 
#'                         annotation.source = "AH",
#'                         genome = "hg19",
#'                         RE = Alu.hg19.demo, 
#'                         ncore = 1,
#'                         verbose = TRUE)
#' }
#' 
#' @export
initREMP <- function(arrayType = c("450k", "EPIC", "Sequencing"), 
                     REtype = c("Alu", "L1", "LTR"),
                     annotation.source = c("AH", "UCSC"), 
                     genome = c("hg19", "hg38"),
                     RE = NULL, 
                     Seq.GR = NULL,
                     ncore = NULL, 
                     BPPARAM = NULL,
                     export = FALSE, 
                     work.dir = tempdir(),
                     verbose = FALSE) {
  ## Initiate running time
  currenT <- Sys.time()

  arrayType <- match.arg(arrayType)
  REtype <- match.arg(REtype)
  annotation.source = match.arg(annotation.source)
  genome = match.arg(genome)
  
  if (!is.null(Seq.GR)) {
    .isGROrStop(Seq.GR)
    names(Seq.GR) <- NULL
  }
  
  message("Start ", REtype, " annotation data initialization (", genome, ")...")
  message("Gene annotation database: ", annotation.source)
  message("Illumina platform: ", arrayType)

  if (is.null(ncore)) ncore <- 1

  ## Setup backend for parallel computing
  be <- getBackend(ncore, BPPARAM, verbose)

  if (arrayType == "450k") {
    if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
      suppressPackageStartupMessages(require("IlluminaHumanMethylation450kanno.ilmn12.hg19"))
      ILMN.GR <- minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      if(genome == "hg38") {
        ILMN.GR <- .liftOver_Hg19toHg38(ILMN.GR, 
                                        "Lifting over Illumina CpG probe location from hg19 to hg38...",
                                        verbose)
      }
    } else stop("Please install missing package: IlluminaHumanMethylation450kanno.ilmn12.hg19")
  } else if (arrayType == "EPIC") {
    if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", quietly = TRUE)) {
      suppressPackageStartupMessages(require("IlluminaHumanMethylationEPICanno.ilm10b2.hg19"))
      ILMN.GR <- minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      if(genome == "hg38") {
        ILMN.GR <- .liftOver_Hg19toHg38(ILMN.GR, 
                                        "Lifting over Illumina CpG probe location from hg19 to hg38...",
                                        verbose)
      }
    } else stop("Please install missing package: IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  } else if (arrayType == "Sequencing") {
    if (!is.null(Seq.GR)) {
      ILMN.GR <- Seq.GR
    } else {
      stop("Seq.GR must be specified if arrayType == 'Sequencing'.")
    }
  }
  
  if (arrayType == "Sequencing") {
    ILMN.GR$Index <- paste0(seqnames(ILMN.GR), ":", start(ILMN.GR))
  } else {
    ILMN.GR <- ILMN.GR[substring(names(ILMN.GR), 1, 2) != "ch"] # remove ch probes
    ILMN.GR$Index <- names(ILMN.GR)
  }

  if (verbose) {
    message(
      "Valid RE-CpG flanking window size: ",
      remp_options(".default.max.flankWindow"), " bp"
    )
  }

  ##################################################
  ### ------ Part I - Resources preparation ------###
  ##################################################

  if (is.null(RE)) {
    ### Get RE annotation database
    RE.hg <- fetchRMSK(REtype = REtype, 
                       annotation.source = annotation.source,
                       genome = genome, 
                       verbose = verbose)
  } else {
    REtype_db <- .guessREtype(RE)
    if(REtype_db != REtype){
      message("Specified REtype (", REtype, ")", " does not match ", REtype_db, " database provided. REtype is set to '", REtype_db, "'.")
      REtype <- REtype_db
    }
    RE.hg <- RE
  }

  ### Get refSeq gene database
  refgene.hg <- fetchRefSeqGene(annotation.source, genome, mainOnly = FALSE, verbose)

  ### Get RE-CpG location database Narrow down RE.hg to RE sequence that
  ### overlaps with CpG sites flanking region For demo data, this will not
  ### change anything.
  ILMN.GR.flank <- .twoWayFlank(ILMN.GR, remp_options(".default.max.flankWindow"))
  RE.hg <- subsetByOverlaps(RE.hg, ILMN.GR.flank, ignore.strand = TRUE)

  ## Narrow down RE.hg to RE sequence that overlaps with gene sequence
  ## RE.hg <- subsetByOverlaps(RE.hg, refgene.hg$main, ignore.strand
  ## = TRUE)

  ## Locate RE-CpG
  RE.CpG <- findRECpG(RE.hg, REtype, genome, be, verbose)

  ## Make RE.hg and RE.CpG dataset consistent with RE
  RE.hg <- RE.hg[base::match(runValue(RE.CpG$Index), runValue(RE.hg$Index))]

  ########################################
  ### ------ Part II - Annotation ------###
  ########################################

  ### Annotate ILMN CpG
  ## Narrow down Illumina CpG Probe using RE.CpG flanking region
  RE.CpG.flanking <- .twoWayFlank(RE.CpG, remp_options(".default.max.flankWindow"))
  ILMN.GR <- subsetByOverlaps(ILMN.GR, RE.CpG.flanking, ignore.strand = TRUE)

  ### ----------------------------------------------------------------------------------------------------
  ### Annotate RE
  RE.refGene <- GRannot(RE.hg, refgene.hg, symbol = FALSE, verbose = verbose)

  ### ----------------------------------------------------------------------------------------------------
  ### RE-CpG covered by ILMN
  RECpG_Platform.hits <- findOverlaps(ILMN.GR, RE.CpG, ignore.strand = TRUE)
  mcols(ILMN.GR)$RE.Index <- Rle(NA)
  ILMN.GR$RE.Index[queryHits(RECpG_Platform.hits)] <- RE.CpG[subjectHits(RECpG_Platform.hits)]$Index

  remparcel <- REMParcel(
    REtype = REtype, 
    genome = genome,
    platform = arrayType,
    RefGene = refgene.hg$main,
    RE = RE.refGene, RECpG = RE.CpG,
    ILMN = ILMN.GR
  )

  message("Done.", .timeTrace(currenT)$t_text)

  if (export) {
    saveParcel(remparcel, work.dir, verbose)
  }

  return(remparcel)
} # End of intREMP


## --------------------------------------
## Internal functions

.guessREtype <- function(RE)
{
  .isGROrStop(RE)
  RE_subfamily <- RE$name
  if(is.null(RE_subfamily)) stop("Please provide the 'name' column specifying the RE subfamily.")
  
  REtype <- c("Alu", "L1", "LTR")
  
  RE_count <- c(sum(grepl("Alu", RE_subfamily)),
                sum(grepl("L1", RE_subfamily)),
                sum(grepl("LTR|ERV", RE_subfamily)))
  
  REtype <- REtype[which.max(RE_count)]
  return(REtype)
}


