#' @title Get BiocParallel back-end
#'
#' @description
#' \code{getBackend} is used to obtain \code{BiocParallel} Back-end to allow 
#' parallel computing.
#'
#' @param ncore Number of cores to run parallel computation. By default max number 
#' of cores available in the machine will be utilized. If \code{ncore = 1}, no 
#' parallel computation is allowed.
#' @param BPPARAM An optional \code{BiocParallelParam} instance determining the 
#' parallel back-end to be used during evaluation. If not specified, default 
#' back-end in the machine will be used.
#' @param verbose Logical parameter. Should the function be verbose?
#' 
#' @return A \code{\link{BiocParallel}} object that can be used for parallel 
#' computing.
#'
#' @examples
#' # Non-parallel mode
#' be <- getBackend(1, verbose = TRUE)
#' be
#' 
#' # parallel mode (2 workers)
#' be <- getBackend(2, verbose = TRUE)
#' be
#' 
#' @export
getBackend <- function(ncore, BPPARAM = NULL, verbose = FALSE) {
  if (ncore > parallel::detectCores()) 
    ncore <- parallel::detectCores()
  if (ncore == 1) {
    backend <- BiocParallel::SerialParam()
    if (verbose) 
      message("You have successfully set non-parallel mode (single worker).")
  } else {
    if (is.null(BPPARAM)) 
      bpparam <- BiocParallel::bpparam() else bpparam <- BPPARAM
      backend_class <- attributes(bpparam)$class
      if (!backend_class %in% names(BiocParallel::registered())) {
        stop("Back-end '", backend_class, "' is not supported on ", 
             Sys.info()["sysname"], "!")
      } else {
        if (backend_class == "SnowParam") 
          backend <- BiocParallel::SnowParam(workers = ncore)
        if (backend_class == "MulticoreParam") 
          backend <- BiocParallel::MulticoreParam(workers = ncore)
        if (verbose) 
          message("You have successfully set parallel mode with ", 
                  BiocParallel::bpworkers(backend), " workers (", backend_class, ").")
      }
  }
  return(backend)
}



#' @title Get RE database from RepeatMasker
#'
#' @description
#' \code{fetchRMSK} is used to obtain specified RE database from RepeatMasker Database 
#' provided by AnnotationHub.
#'
#' @param ah An \code{\link{AnnotationHub}} object. Use \code{AnnotationHub()} to 
#' retrive information about all records in the hub.
#' @param REtype Type of RE. Currently \code{"Alu"} and \code{"L1"} are supported.
#' @param verbose Logical parameter. Should the function be verbose?
#' 
#' @return A \code{\link{GRanges}} object containing RE database. 'name' column 
#' indicates the RE subfamily; 'score' 
#' column indicates the SW score; 'Index' is an internal index for RE to facilitate data 
#' referral, which is meaningless for external use.
#'
#' @examples
#' ah <- AnnotationHub::AnnotationHub()
#' L1 <- fetchRMSK(ah, 'L1', verbose = TRUE)
#' L1
#' 
#' @export
fetchRMSK <- function(ah, REtype, verbose = FALSE) {
  if (verbose) {
    message("Loading ", REtype, " annotation data from RepeatMasker (hg19) database from AnnotationHub: ", 
            remp_options(".default.AH.repeatmasker.hg19"))
  }
  # HG19 rtracklayer://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk
  rmsk <- suppressWarnings(suppressMessages(ah[[remp_options(".default.AH.repeatmasker.hg19")]]))  
  
  if (REtype == "Alu") 
    REFamily <- remp_options(".default.AluFamily")
  if (REtype == "L1") 
    REFamily <- remp_options(".default.L1Family")
  
  RE <- rmsk[rmsk$name %in% REFamily, ]
  RE <- RE[seqnames(RE) %in% paste0("chr", c(1:22, "X", "Y"))]  # chr1 - 22, chrX, chrY
  
  mcols(RE)$Index <- Rle(paste(REtype, formatC(seq_len(length(RE)), width = 7, 
                                               flag = "0"), sep = "_"))  # Add internal index, meaningless for external database
  return(RE)
}



#' @title Get RefSeq gene database
#'
#' @description
#' \code{fetchRefSeqGene} is used to obtain refSeq gene database provided by AnnotationHub.
#'
#' @param ah An \code{AnnotationHub} object. Use \code{AnnotationHub()} to retrive information 
#' about all records in the hub.
#' @param mainOnly Logical parameter. See details.
#' @param verbose Logical parameter. Should the function be verbose?
#'
#' @details 
#' When \code{mainOnly = FALSE}, only the transcript location information will be returned, 
#' Otherwise, a \code{\link{GRangesList}} object containing gene regions 
#' information will be added. Gene regions include: 2000 base pair upstream of the transcript 
#' start site (\code{$tss})), 5'UTR (\code{$fiveUTR})), coding sequence (\code{$cds})), 
#' exon (\code{$exon})), and 3'UTR (\code{$threeUTR})). The \code{index} column is an internal 
#' index that is used to facilitate data referral, which is meaningless for external use.
#' 
#' @return A single \code{\link{GRanges}} (for main refgene data) object or a list incorporating 
#' both \code{GRanges} object (for main refgene data)  and \code{\link{GRangesList}} object 
#' (for gene regions data).
#' 
#' @examples
#' ah <- AnnotationHub::AnnotationHub()
#' refGene <- fetchRefSeqGene(ah, mainOnly = TRUE, verbose = TRUE)
#' refGene
#' 
#' @export
fetchRefSeqGene <- function(ah, mainOnly = FALSE, verbose = FALSE) {
  if (verbose) {
    message("Loading refseq gene (hg19) database from AnnotationHub: ", 
            remp_options(".default.AH.refgene.hg19"))
    
  }
  
  refgene <- suppressWarnings(suppressMessages(ah[[remp_options(".default.AH.refgene.hg19")]]))  
  refgene <- refgene[as.character(seqnames(refgene)) %in% 
                       paste0("chr", c(1:22, "X", "Y"))]  ## chr1 - 22, chrX, chrY
  refgene$type <- substring(refgene$name, 1, 2)  ## Protein coding gene (NM) or noncoding RNA gene (NR)
  refgene$index <- as.character(seq_len(length(refgene)))  ## This is internal index, meaningless for external database
  
  ## Master database with detailed gene information
  refgene_main <- GRanges(seqnames = seqnames(refgene), ranges = ranges(refgene), 
                          strand = strand(refgene), name = mcols(refgene)$name, type = mcols(refgene)$type, 
                          index = mcols(refgene)$index)
  
  refseq2EG <- unlist(BiocGenerics::mget(refgene_main$name, org.Hs.eg.db::org.Hs.egREFSEQ2EG, 
                                         ifnotfound = NA))
  if (verbose) {
    message("Note: There are ", sum(is.na(refseq2EG)), " refseq ID cannot be mapped to Entrez ID, which could be obsolete.")
  }
  # obsoleteID <- refgene_main$name[is.na(refseq2EG)]
  refseq2EG <- na.omit(refseq2EG)
  
  EG2Symbol <- unlist(BiocGenerics::mget(refseq2EG, org.Hs.eg.db::org.Hs.egSYMBOL))  ## If EG -> symbol not found, this is a fatel error, should stop here
  
  refseq2EG_Symbol.DF <- DataFrame(refseqID = names(refseq2EG), EntrezGeneID = refseq2EG, 
                                   GeneSymbol = EG2Symbol)
  
  refgene_main$EntrezGene <- refseq2EG_Symbol.DF$EntrezGeneID[match(refgene_main$name, 
                                                                    refseq2EG_Symbol.DF$refseqID)]
  refgene_main$GeneSymbol <- refseq2EG_Symbol.DF$GeneSymbol[match(refgene_main$name, 
                                                                  refseq2EG_Symbol.DF$refseqID)]
  
  if (mainOnly)
    return(refgene_main)
  
  ## TSS2000
  refgene_tss <- promoters(refgene, upstream = remp_options(".default.TSS.upstream"), 
                           downstream = remp_options(".default.TSS.downstream"))[, "index"]
  
  ## CDS
  refgene_cds <- GRanges(seqnames = seqnames(refgene), ranges = refgene$thick, 
                         strand = strand(refgene), index = mcols(refgene)$index)
  
  ## EXON
  block_shift <- shift(refgene$blocks, start(refgene) - 1)  # minus one
  refgene_exon <- makeGRangesListFromFeatureFragments(seqnames = seqnames(refgene_main), 
                                                      fragmentStarts = start(block_shift), fragmentEnds = end(block_shift), 
                                                      strand = strand(refgene_main))
  names(refgene_exon) <- refgene_main$index
  refgene_exon <- unlist(refgene_exon)
  refgene_exon$index <- as.integer(names(refgene_exon))
  names(refgene_exon) <- NULL
  
  ## UTR (5'UTR/3'UTR concept only applied for protein-coding gene)
  refgene_NM <- refgene[refgene$type == "NM", c("thick", "index")]  
  refgene_f <- refgene_NM[strand(refgene_NM) == "+"]
  refgene_r <- refgene_NM[strand(refgene_NM) == "-"]
  
  refgene_5UTR_f <- GRanges(seqnames = seqnames(refgene_f), 
                            ranges = IRanges(start = start(refgene_f), 
                                             end = start(refgene_f$thick)), 
                            strand = strand(refgene_f), 
                            index = refgene_f$index)
  
  refgene_5UTR_r <- GRanges(seqnames = seqnames(refgene_r), 
                            ranges = IRanges(start = end(refgene_r$thick), 
                                             end = end(refgene_r)), 
                            strand = strand(refgene_r), 
                            index = refgene_r$index)
  
  refgene_3UTR_f <- GRanges(seqnames = seqnames(refgene_f), 
                            ranges = IRanges(start = end(refgene_f$thick), 
                                             end = end(refgene_f)), 
                            strand = strand(refgene_f), 
                            index = refgene_f$index)
  
  refgene_3UTR_r <- GRanges(seqnames = seqnames(refgene_r), 
                            ranges = IRanges(start = start(refgene_r), 
                                             end = start(refgene_r$thick)), 
                            strand = strand(refgene_r), 
                            index = refgene_r$index)
  
  refgene_5UTR <- c(refgene_5UTR_f, refgene_5UTR_r)
  refgene_3UTR <- c(refgene_3UTR_f, refgene_3UTR_r)
  
  return(list(main = refgene_main, 
              regions = GRangesList(tss = refgene_tss, 
                                    cds = refgene_cds, 
                                    exon = refgene_exon, 
                                    fiveUTR = refgene_5UTR, 
                                    threeUTR = refgene_3UTR))
  )
  
}



#' @title Find RE-CpG genomic location given RE ranges information
#'
#' @description
#' \code{findRECpG} is used to obtain RE-CpG genomic location data.
#'
#' @param RE.hg19 an \code{GRanges} object of RE genomic location database. This 
#' can be obtained by \code{\link{fetchRMSK}}.
#' @param REtype Type of RE. Currently \code{"Alu"} and \code{"L1"} are supported.
#' @param be A \code{\link{BiocParallel}} object containing back-end information that is 
#' ready for paralle computing. This can be obtained by \code{\link{getBackend}}.
#' @param verbose logical parameter. Should the function be verbose?
#' 
#' @details 
#' CpG site is defined as 5'-C-p-G-3'. It is reasonable to assume that the methylation 
#' status across all CpG/CpG dyads are concordant. Maintenance methyltransferase exhibits 
#' a preference for hemimethylated CpG/CpG dyads (methylated on one strand only). 
#' As a result, methyaltion status of CpG sites in both forward and reverse strands are usually consistent. 
#' Therefore, to accommodate the cytosine loci in both strands, the returned genomic 
#' ranges cover the 'CG' sequence with width of 2. The 'stand' information indicates the strand of the RE.
#' 
#' @return A \code{\link{GRanges}} object containing identified RE-CpG genomic 
#' location data.
#'
#' @examples
#' data(Alu.demo)
#' be <- getBackend(1, verbose = TRUE)
#' RE.CpG <- findRECpG(Alu.demo, 'Alu', be, verbose = TRUE) 
#' RE.CpG
#' 
#' @export
findRECpG <- function(RE.hg19, REtype = c("Alu", "L1"), be = NULL, verbose = FALSE) {
  
  if(is.null(be)) be <- getBackend(1, verbose)
  REtype <- match.arg(REtype)
  
  ## Flank RE by 1 base pair to cover the cases where CG is located in
  ## RE's head or tail
  RE.hg19 <- .twoWayFlank(RE.hg19, 1)
  
  if (verbose) 
    message("    Getting sequence of the total ", length(RE.hg19), 
            " ", REtype, " ...")
  SEQ.RE <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, 
                             GRanges(seqnames(RE.hg19), 
                                     IRanges(start = start(RE.hg19), 
                                             end = end(RE.hg19))))  # ignore strand
  
  if (verbose) 
    message("    Identifying CpG sites in ", REtype, " sequence ...")
  
  if (BiocParallel::bpworkers(be) > 1) {
    bpstart(be)
    .bploadLibraryQuiet("Biostrings", be)
    RE.CpG <- bpvec(SEQ.RE, .vRECpGPos, CpG = DNAString("CG"), BPPARAM = be)
    bpstop(be)
  } else {
    RE.CpG <- .vRECpGPos(SEQ.RE, CpG = DNAString("CG"))
  }
  
  noCpGind <- which(sapply(RE.CpG, length) == 0)  ## Remove sequence without CpG sites
  
  if (length(noCpGind) == 0) {
    if (verbose) 
      message("    All ", REtype, " sequences contain CpGs.")
  } else {
    if (verbose) 
      message("    There are ", length(noCpGind), " ", REtype, " sequences contain no CpG site.")
    RE.hg19 <- RE.hg19[-noCpGind]
    RE.CpG <- RE.CpG[-noCpGind]
  }
  
  ## Extend 'C' to 'CG' AND Shift to the real genomic location
  RE.CpG.IRList <- IRangesList(start = RE.CpG - 1, end = RE.CpG)
  block_shift <- shift(RE.CpG.IRList, start(RE.hg19))
  
  RE.CpG <- makeGRangesListFromFeatureFragments(seqnames = seqnames(RE.hg19), 
                                                fragmentStarts = start(block_shift), 
                                                fragmentEnds = end(block_shift), 
                                                strand = strand(RE.hg19))
  elementLen <- elementNROWS(RE.CpG)
  
  RE.CpG <- unlist(RE.CpG)
  RE.CpG$Index <- Rle(runValue(RE.hg19$Index), elementLen)
  RE.CpG <- RE.CpG[order(RE.CpG$Index)]
  
  if (verbose) 
    message("    Identified ", length(RE.CpG), " CpG/CpG dyads in ", 
            nrun(RE.CpG$Index), " ", REtype)
  
  return(RE.CpG)
}




#' @title Annotate genomic ranges data with gene region information.
#'
#' @description
#' \code{GRannot} is used to annotate a \code{GRanges} dataset with gene region 
#' information using refseq gene database
#'
#' @param object.GR An \code{\link{GRanges}} object of a genomic location database.
#' @param refgene.hg19 A complete refGene annotation database returned by 
#' \code{\link{fetchRefSeqGene}} (with parameter \code{mainOnly = FALSE}).
#' @param verbose Logical parameter. Should the function be verbose?
#' 
#' @details
#' The annotated gene region information includes: protein coding gene (InNM), 
#' noncoding RNA gene (InNR), 2000 base pair upstream of the transcript start site (InTSS),
#' 5'UTR (In5UTR), coding sequence (InCDS), exon (InExon), and 3'UTR (In3UTR). The intergenic 
#' and intron regions can then be represented by the combination of these region data. 
#' The number shown in these columns represent the row number or 'index' column in the 
#' main refgene database obtained by \code{\link{fetchRefSeqGene}}.
#' 
#' @return A \code{\link{GRanges}} or a \code{\link{GRangesList}} object containing refSeq 
#' Gene database.
#'
#' @examples
#' data(Alu.demo)
#' ah <- AnnotationHub::AnnotationHub()
#' refgene.hg19 <- fetchRefSeqGene(ah, verbose = TRUE)
#' Alu.demo.refGene <- GRannot(Alu.demo, refgene.hg19, verbose = TRUE)
#' Alu.demo.refGene
#' 
#' @export
GRannot <- function(object.GR, refgene.hg19, verbose = FALSE) {
  
  ## Check if the object.GR contains a metadata column called "Index"
  if(!"Index" %in% colnames(mcols(object.GR)))
    object.GR$Index <- as.character(seq_len(length(object.GR)))
  
  if(any(duplicated(as.character(object.GR$Index))))
    stop("The 'Index' column provided in the GRanges object must be unique.")
  
  ####### InNM
  NM.GR <-  refgene.hg19$main[refgene.hg19$main$type == "NM"]
  NM.GR <- .oneWayFlank(NM.GR, 2000, start = TRUE)
  object.GR <- .addannot(object.GR, 
                         NM.GR, 
                         "InNM", 
                         verbose)
  
  ####### InNR
  NR.GR <-  refgene.hg19$main[refgene.hg19$main$type == "NR"]
  NR.GR <- .oneWayFlank(NR.GR, 2000, start = TRUE)
  object.GR <- .addannot(object.GR, 
                         NR.GR, 
                         "InNR", 
                         verbose)
  
  ####### InTSS2000
  object.GR <- .addannot(object.GR, 
                         refgene.hg19$regions$tss, 
                         "InTSS", 
                         verbose)
  
  ####### In5UTR
  object.GR <- .addannot(object.GR, 
                         refgene.hg19$regions$fiveUTR, 
                         "In5UTR", 
                         verbose)
  
  ####### InCDS
  object.GR <- .addannot(object.GR, 
                         refgene.hg19$regions$cds, 
                         "InCDS", 
                         verbose)
  
  ####### InExon
  object.GR <- .addannot(object.GR, 
                         refgene.hg19$regions$exon, 
                         "InExon", 
                         verbose)
  
  ####### In3UTR
  object.GR <- .addannot(object.GR, 
                         refgene.hg19$regions$threeUTR, 
                         "In3UTR", 
                         verbose)
  
  return(object.GR)
}

##--------------------------------------
## Internal functions

# Used by findRECpG (vetorized function)
.vRECpGPos <- function(seq, CpG) {
  start(vmatchPattern(CpG, seq))
}

# Used by addannot
.clps <- function(x) paste0(x, collapse = "|")

# Used by GRannot
.addannot <- function(object.GR, annot.GR, annotName, verbose) {
  if (verbose) 
    message("    ", annotName, " ... ", appendLF = FALSE)
  
  object_annot_Hit <- findOverlaps(object.GR, annot.GR, ignore.strand = TRUE)
  object_raw <- DataFrame(InRegion = annot.GR[subjectHits(object_annot_Hit)]$index)
  object_agg <- aggregate(object_raw, list(Index = as.character(object.GR$Index[queryHits(object_annot_Hit)])), 
                          .clps)
  
  InRegion <- rep(NA, length(object.GR))
  InRegion[match(object_agg$Index, object.GR$Index)] <- object_agg$InRegion
  
  object.GR$newRegion <- InRegion
  mcols(object.GR) <- .changeColNames(mcols(object.GR), "newRegion", 
                                      annotName)
  
  if (verbose) 
    message("  ", nrow(object_agg), " found!")
  
  return(object.GR)
}