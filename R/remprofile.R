#' @title Extract DNA methylation data profiled in RE
#'
#' @description
#' \code{remprofile} is used to extract profiled methylation of CpG sites in RE.
#'
#' @param methyDat A \code{\link{RatioSet}}, \code{\link{GenomicRatioSet}}, \code{\link{DataFrame}},
#' \code{data.table}, \code{data.frame}, or \code{matrix} of Illumina BeadChip methylation data
#' (450k or EPIC array) or Illumina methylation percentage estimates by sequencing.
#' @param genome Character parameter. Specify the build of human genome. Can be either \code{"hg19"} or 
#' \code{"hg38"}. For 450k/EPIC array, \code{"hg19"} is used more often while specifying \code{"hg38"} will 
#' lift over the Illumina CpG probe location to build \code{"hg38"}. For sequencing data, please make
#' sure the specified genome build is consistent with the actual genome build of \code{Seq.GR}.
#' @param REtype Type of RE. Currently \code{"Alu"}, \code{"L1"}, and \code{"LTR"} are supported.
#' @param annotation.source Character parameter. Specify the source of annotation databases, including
#' the RefSeq Gene annotation database and RepeatMasker annotation database. If \code{"AH"}, the database 
#' will be obtained from the AnnotationHub package. If \code{"UCSC"}, the database will be downloaded 
#' from the UCSC website http://hgdownload.cse.ucsc.edu/goldenpath. The corresponding build (\code{"hg19"} or 
#' \code{"hg38"}) will be specified in the parameter \code{genome}.
#' @param Seq.GR A \code{\link{GRanges}} object containing genomic locations of the CpGs profiled by sequencing
#' platforms. This parameter should not be \code{NULL} if the input methylation data \code{methyDat} are
#' obtained by sequencing. Note that the genomic location can be in either hg19 or hg38 build. The user
#' should make sure the parameter \code{genome} is correctly specified. 
#' @param RE A \code{\link{GRanges}} object containing user-specified RE genomic location information.
#' If \code{NULL}, the function will retrive RepeatMasker RE database from \code{\link{AnnotationHub}}
#' (build hg19) or download the database from UCSC website (build hg19/hg38).
#' @param impute Parameter used by \code{\link{grooMethy}}. If \code{TRUE}, K-Nearest Neighbouring
#' imputation will be applied to fill the missing values. Default = \code{FALSE}.
#' @param imputebyrow Parameter used by \code{\link{grooMethy}}. If \code{TRUE}, missing values will
#' be imputed using similar values in row (i.e., across samples); if \code{FALSE}, missing values
#' will be imputed using similar values in column (i.e., across CpGs). Default is \code{TRUE}.
#' @param verbose Logical parameter. Should the function be verbose?
#'
#' @return A \code{\link{REMProduct}} object containing profiled RE methylation results.
#'
#' @examples
#' data(Alu.hg19.demo)
#' if (!exists("GM12878_450k")) GM12878_450k <- getGM12878("450k")
#' remprofile.res <- remprofile(GM12878_450k,
#'                              REtype = "Alu",
#'                              annotation.source = "AH",
#'                              genome = "hg19",
#'                              RE = Alu.hg19.demo,
#'                              verbose = TRUE)
#' details(remprofile.res)
#' rempB(remprofile.res) # Methylation data (beta value)
#' 
#' remprofile.res <- rempAggregate(remprofile.res)
#' details(remprofile.res)
#' rempB(remprofile.res) # Methylation data (beta value)
#' 
#' @export
remprofile <- function(methyDat, 
                       REtype = c("Alu", "L1", "LTR"), 
                       annotation.source = c("AH", "UCSC"), 
                       genome = c("hg19", "hg38"), 
                       Seq.GR = NULL,
                       RE = NULL, 
                       impute = FALSE, 
                       imputebyrow = TRUE, 
                       verbose = FALSE) {
  if (is.null(methyDat)) stop("Methylation dataset (methyDat) is missing.")
  if (!is.null(Seq.GR)) {
    .isGROrStop(Seq.GR)
    names(Seq.GR) <- NULL
  }
  
  genome <- match.arg(genome)
  REtype <- match.arg(REtype)
  annotation.source = match.arg(annotation.source)
  
  ## Groom methylation data
  methyDat <- grooMethy(methyDat, Seq.GR, impute, verbose = verbose)
  arrayType <- gsub("IlluminaHumanMethylation", "", methyDat@annotation["array"])
  
  methyDat <- minfi::getM(methyDat)
  probeNames <- rownames(methyDat)
  
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

  probeNames <- intersect(probeNames, ILMN.GR$Index)

  methyDat <- methyDat[probeNames, , drop = FALSE]
  ILMN.GR <- ILMN.GR[base::match(probeNames, ILMN.GR$Index)]
  # identical(names(ILMN.GR), rownames(methyDat))

  if (is.null(RE)) {
    ### Get RE annotation database
    RE.hg <- fetchRMSK(REtype = REtype, 
                       annotation.source = annotation.source,
                       genome = genome, 
                       verbose = verbose)
  } else {
    .isGROrStop(RE)
    RE.hg <- RE
  }

  ### RE-CpG covered by ILMN
  RECpG_Platform.hits <- findOverlaps(RE.hg, ILMN.GR, ignore.strand = TRUE)

  # Update RE ranges
  RE.hg <- RE.hg[queryHits(RECpG_Platform.hits)]

  # Update CpG ranges
  mcols(ILMN.GR)$RE.Index <- Rle(NA)
  ILMN.GR$RE.Index[subjectHits(RECpG_Platform.hits)] <- RE.hg$Index
  RE_CpG_ILMN <- ILMN.GR[!is.na(ILMN.GR$RE.Index)]
  RE_CpG_ILMN <- RE_CpG_ILMN[order(RE_CpG_ILMN$RE.Index)]

  cpgRanges <- RE_CpG_ILMN

  RE_annotation <- unique(RE.hg)
  refgene.hg <- fetchRefSeqGene(annotation.source = annotation.source, genome = genome, mainOnly = FALSE, verbose)
  refgene_main <- refgene.hg$main

  RE_annotation <- GRannot(RE_annotation, refgene.hg, symbol = FALSE, verbose = verbose)

  RE_annotation_name <- colnames(mcols(RE_annotation))
  regionCode <- mcols(RE_annotation)[remp_options(".default.genomicRegionColNames")]
  RE_annotation <- RE_annotation[, RE_annotation_name[!RE_annotation_name %in%
    remp_options(".default.genomicRegionColNames")]]

  samplenames <- colnames(methyDat)
  sampleN <- length(samplenames)

  sampleinfo <- DataFrame(matrix(NA, nrow = sampleN, ncol = 1))
  colnames(sampleinfo) <- "Not_Applicable"
  rownames(sampleinfo) <- samplenames

  Profiled_RECpG_M <- methyDat[cpgRanges$Index, , drop = FALSE]

  ## RE coverage
  RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN,
    REtype,
    indent = "    ", verbose
  )

  # Gene coverage
  GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main,
    REtype,
    indent = "    ", verbose
  )

  remproduct <- REMProduct(
    REtype = REtype, genome = genome,
    platform = arrayType, win = "N/A",
    predictModel = "Profiled", QCModel = "N/A",
    rempM = Profiled_RECpG_M, rempQC = NULL,
    cpgRanges = cpgRanges, sampleInfo = sampleinfo,
    REannotation = RE_annotation,
    RECpG = ILMN.GR,
    regionCode = regionCode,
    refGene = refgene_main,
    varImp = DataFrame(),
    REStats = RE_COVERAGE, GeneStats = GENE_COVERAGE,
    Seed = NA
  )
  return(remproduct)
}
