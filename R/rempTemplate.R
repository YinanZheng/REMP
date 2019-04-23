#' @title Prepare data template for REMP
#'
#' @description
#' \code{rempTemplate} is used to build a set of data templates for prediction. The data templates include
#' RE-CpGs and their methylation data for model training, neighboring CpGs of RE-CpGs and their
#' methylation data for model prediction, and other necessary information about the prediction.
#' This function is useful when one needs to experiment different tunning parameters so that these pre-built
#' data templates can be re-used.
#'
#' @param methyDat A \code{\link{RatioSet}}, \code{\link{GenomicRatioSet}}, \code{\link{DataFrame}},
#' \code{data.table}, \code{data.frame}, or \code{matrix} of Illumina BeadChip methylation data
#' (450k or EPIC array).
#' @param remparcel An \code{\link{REMParcel}} object containing necessary data to carry out the
#' prediction. If \code{NULL}, the function will search the \code{.rds} data file in \code{work.dir}
#' exported by \code{\link{initREMP}} (with \code{export = TRUE}) or \code{\link{saveParcel}}.
#' @param win An integer specifying window size to confine the upstream and downstream flanking
#' region centered on the predicted CpG in RE for prediction. Default = \code{1000}.
#' @param verbose Logical parameter. Should the function be verbose?
#'
#' @return A \code{template} object containing a \code{\link{GRanges}} object of neifhboring CpGs of RE-CpGs
#' to be predicted (\code{$NBCpG_GR}) and their methylation dataset matrix (\code{$NBCpG_methyDat});
#' a \code{\link{GRanges}} object of RE-CpGs for model training (\code{$RECpG_GR}) and their methylation
#' dataset matrix (\code{$RECpG_methyDat}); \code{\link{GRanges}} objects of RefSeq Gene database (\code{$refgene})
#' and RE (\code{$RE}) that are extracted from the \code{parcel} input; a string of RE type (\code{$REtype})
#' and a string of methylation platform (\code{$arrayType}). Note: the subset operator \code{[]} is supported.
#'
#' @examples
#' if (!exists("GM12878_450k")) GM12878_450k <- getGM12878("450k")
#' if (!exists("remparcel")) {
#'   data(Alu.demo)
#'   remparcel <- initREMP(arrayType = "450k", REtype = "Alu", RE = Alu.demo, ncore = 1)
#' }
#' 
#' template <- rempTemplate(GM12878_450k, remparcel, win = 1000)
#' template
#' 
#' ## To make a subset
#' template[1]
#' @export
rempTemplate <- function(methyDat = NULL, remparcel = NULL, win = 1000, verbose = FALSE) {
  if (is.null(methyDat)) stop("Methylation dataset is missing.")

  if (is.null(remparcel)) stop("REMParcel object is missing.")
  .isREMParcelOrStop(remparcel)

  if (win > remp_options(".default.max.flankWindow")) {
    stop(
      "Flanking window size cannot be greater than ", remp_options(".default.max.flankWindow"),
      ". Please see ?remp_options for details."
    )
  }

  t <- Sys.time()

  ## Groom methylation data
  methyDat <- grooMethy(methyDat, verbose = verbose)
  methyDat <- minfi::getM(methyDat)

  ## Guess array type
  arrayType <- remparcel@REMParcelInfo$platform
  REtype <- remparcel@REMParcelInfo$REtype
  RE_refGene.original <- getRE(remparcel)
  RE_CpG <- getRECpG(remparcel)
  ILMN <- getILMN(remparcel)
  RE_CpG_ILMN <- getILMN(remparcel, REonly = TRUE)

  ## Make indicator of genomic regions
  RE_refGene <- .toIndicator(RE_refGene.original)

  ## RE.CpG : Add CpG density and rename
  ## identical(runValue(RE_refGene$Index), runValue(RE_CpG$Index))
  RE_refGene$N <- runLength(RE_CpG$Index)
  RE_refGene$Length <- width(RE_refGene)
  RE_refGene$CpG.density <- RE_refGene$N / RE_refGene$Length # density of CpG within RE sequence

  merge_RE_CpG_meta <- mcols(RE_refGene)[match(
    as.character(RE_CpG$Index),
    as.character(RE_refGene$Index)
  ), ]
  merge_RE_CpG_meta$CpG.ID <- as.character(granges(RE_CpG))
  mcols(RE_CpG) <- setNames(merge_RE_CpG_meta, paste0("RE.", colnames(merge_RE_CpG_meta)))
  # all(as.character(granges(RE_CpG)) == RE_CpG$RE.CpG.ID)

  ## Find valid CpG sites with methylation profiled
  commonCpGIndex <- intersect(rownames(methyDat), ILMN$Index)

  ## Trim down methyDat:
  methyDat <- methyDat[match(commonCpGIndex, rownames(methyDat)), , drop = FALSE]

  ## ILMN : create pointer and rename. Note: Different Methy data will
  ## result in different size of this data!
  ILMN <- ILMN[match(commonCpGIndex, ILMN$Index), ]
  mcols(ILMN) <- setNames(mcols(ILMN), paste0("ILMN.", colnames(mcols(ILMN))))
  ILMN$Methy.ptr <- seq_len(length(ILMN))
  # identical(ILMN$ILMN.Index, rownames(methyDat))

  ## RE_CpG_ILMN : Compute total number of RE covered by ILMN. Note:
  ## Different Methy data will result in different size of this data.
  RE_CpG_ILMN <- RE_CpG_ILMN[RE_CpG_ILMN$Index %in% rownames(methyDat), ]
  RE_CpG_ILMN <- unique(RE_CpG_ILMN)

  RE_CpG_ILMN_DATA <- methyDat[match(RE_CpG_ILMN$Index, rownames(methyDat)), , drop = FALSE]

  ############################################################# Find RE_CpG's neighboring CpG
  message("Processing ", REtype, " with window +/- ", win, " base pair ...")

  RE_CpG_flanking <- .twoWayFlank(granges(RE_CpG), win)
  HITS <- findOverlaps(RE_CpG_flanking, ILMN, ignore.strand = TRUE)

  ## Part 1: RE-CpG
  RE_NeibCpG <- RE_CpG[queryHits(HITS), ]
  ## Total RE that can be predicted (contains neighboring CpGs within given window)

  ## Part 2: Neighboring ILMN CpG
  RE_NeibCpG_ILMN <- ILMN[subjectHits(HITS), ]
  RE_NeibCpG_ILMN <- DataFrame(
    RE_NeibCpG_ILMN.GR = granges(RE_NeibCpG_ILMN),
    mcols(RE_NeibCpG_ILMN)
  )
  mcols(RE_NeibCpG) <- DataFrame(mcols(RE_NeibCpG), RE_NeibCpG_ILMN)
  # identical(as.character(granges(RE_NeibCpG)), RE_NeibCpG$RE.CpG.ID)

  ## Remove singleton (RE-CpGs with only one neighboring ILMN CpG)
  RE_NeibCpG$distance <- abs(start(RE_NeibCpG$RE_NeibCpG_ILMN.GR) - start(RE_NeibCpG))
  RE_NeibCpG <- RE_NeibCpG[RE_NeibCpG$distance > 1, ]

  ## Add core predictors.
  RE_NeibCpG <- RE_NeibCpG[order(RE_NeibCpG$RE.CpG.ID, RE_NeibCpG$distance), ]
  RE_NeibCpG.DF <- mcols(RE_NeibCpG)[, c("RE.CpG.ID", "Methy.ptr", "distance")]

  distance_agg <- .aggregateNeib(
    distance ~ RE.CpG.ID, RE_NeibCpG.DF,
    function(x) c(length(x), mean(x), sd(x), x[1], x[2]),
    c(
      "RE.CpG.ID", "N.nbr", "distance.mean",
      "distance.std", "distance.min", "distance.min2"
    )
  )
  Methy_ptr_agg <- .aggregateNeib(
    Methy.ptr ~ RE.CpG.ID, RE_NeibCpG.DF,
    function(x) c(x[1], x[2]),
    c("RE.CpG.ID", "Methy.ptr.min", "Methy.ptr.min2")
  )

  Methy_ptr_distance_agg <- merge(distance_agg, Methy_ptr_agg, by = "RE.CpG.ID")

  RE_NeibCpG_meta <- mcols(RE_NeibCpG)
  mcols(RE_NeibCpG) <- cbind(
    RE_NeibCpG_meta,
    Methy_ptr_distance_agg[match(
      RE_NeibCpG_meta$RE.CpG.ID,
      Methy_ptr_distance_agg$RE.CpG.ID
    ), -1]
  )
  RE_NeibCpG <- RE_NeibCpG[order(RE_NeibCpG$RE.Index)]
  RE_NeibCpG <- RE_NeibCpG[RE_NeibCpG$N.nbr > 1, ] ## Remove singleton

  refgene_main <- getRefGene(remparcel)
  RE_refGene.original <- getRE(remparcel)

  res <- list(
    NBCpG_methyDat = DataFrame(methyDat), RECpG_methyDat = DataFrame(RE_CpG_ILMN_DATA),
    NBCpG_GR = RE_NeibCpG, RECpG_GR = RE_CpG_ILMN,
    refgene = getRefGene(remparcel),
    RE = getRE(remparcel),
    REtype = REtype,
    arrayType = arrayType
  )

  class(res) <- "template"

  message("Done.", .timeTrace(t)$t_text)

  return(res)
}

#' @export
"[.template" <- function(x, i) {
  x$NBCpG_methyDat <- x$NBCpG_methyDat[, i, drop = FALSE]
  x$RECpG_methyDat <- x$RECpG_methyDat[, i, drop = FALSE]
  x
}
