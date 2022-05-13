#' @title Get BiocParallel back-end
#'
#' @description
#' \code{getBackend} is used to obtain \code{BiocParallel} Back-end to allow
#' parallel computing.
#'
#' @param ncore Number of cores used for parallel computing. By default max number
#' of cores available in the machine will be utilized. If \code{ncore = 1}, no
#' parallel computing is allowed.
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
#' be <- getBackend(ncore = 1, verbose = TRUE)
#' be
#' 
#' # parallel mode (2 workers)
#' be <- getBackend(ncore = 2, verbose = TRUE)
#' be
#' 
#' @export
getBackend <- function(ncore, 
                       BPPARAM = NULL, 
                       verbose = FALSE) {
  maxCore <- parallel::detectCores()
  if (ncore > maxCore) {
    if (!is.na(maxCore)) {
      ncore <- maxCore
    } else {
      
    }
  }

  if (ncore > 16) message("Note: Requesting more than 16 cores may not benifit much.")

  if (ncore == 1) {
    backend <- BiocParallel::SerialParam()
    if (verbose) {
      message("You have successfully set non-parallel mode (single worker).")
    }
  } else {
    if (is.null(BPPARAM)) bpparam <- BiocParallel::bpparam() else bpparam <- BPPARAM
    backend_class <- attributes(bpparam)$class
    if (!backend_class %in% names(BiocParallel::registered())) {
      stop(
        "Back-end '", backend_class, "' is not supported on ",
        Sys.info()["sysname"], "!"
      )
    } else {
      if ("MulticoreParam" %in% names(BiocParallel::registered())) {
        backend <- BiocParallel::MulticoreParam(workers = ncore)
      } else {
        backend <- BiocParallel::SnowParam(workers = ncore)
      }
      if (verbose) {
        message(
          "You have successfully set parallel mode with ",
          BiocParallel::bpnworkers(backend), " workers (", backend_class, ")."
        )
      }
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
#' @param REtype Type of RE. Currently \code{"Alu"}, \code{"L1"}, and \code{"ERV"} are supported.
#' @param annotation.source Character parameter. Specify the source of annotation databases, including
#' the RefSeq Gene annotation database and RepeatMasker annotation database. If \code{"AH"}, the database 
#' will be obtained from the AnnotationHub package. If \code{"UCSC"}, the database will be downloaded 
#' from the UCSC website http://hgdownload.cse.ucsc.edu/goldenpath. The corresponding build (\code{"hg19"} or 
#' \code{"hg38"}) will be specified in the parameter \code{genome}.
#' @param genome Character parameter. Specify the build of human genome. Can be either \code{"hg19"} or
#' \code{"hg38"}.
#' @param verbose Logical parameter. Should the function be verbose?
#'
#' @return A \code{\link{GRanges}} object containing RE database. 'repName' column
#' indicates the RE name; 'swScore' column indicates the SW score; 'Index' is an
#' internal index for RE to facilitate data referral, which is meaningless for external use.
#'
#' @examples
#' L1 <- fetchRMSK(REtype = "L1", 
#'                 annotation.source = "AH",
#'                 genome = "hg19", 
#'                 verbose = TRUE)
#' L1
#' @export
fetchRMSK <- function(REtype = c("Alu", "L1", "ERV"), 
                      annotation.source = c("AH", "UCSC"), 
                      genome = c("hg19", "hg38"), 
                      verbose = FALSE) {
  REtype = match.arg(REtype)
  annotation.source = match.arg(annotation.source)
  genome = match.arg(genome)
  
  if(genome == "hg19") {
    if(annotation.source == "UCSC") {
      message("UCSC RepeatMasker annotation database is not recommended. Switching to AnnotationHub instead...")
      annotation.source = "AH"
    }
    if(annotation.source == "AH") {
      if (verbose) message(
        "Loading ", REtype, " annotation data from RepeatMasker (hg19) database from AnnotationHub: ",
        remp_options(".default.AH.repeatmasker.hg19")
      )
      ah <- .initiateAH()
      if(is.null(ah)) {
        warning("AnnotationHub is currently not accessible, Please try again later.")
      } else {
        rmsk <- suppressMessages(ah[[remp_options(".default.AH.repeatmasker.hg19")]])
      }
    }
  }
  
  if(genome == "hg38") {
    if(annotation.source == "UCSC") {
      message("UCSC RepeatMasker annotation database is not recommended. Switching to AnnotationHub instead...")
      annotation.source = "AH"
    }
    if(annotation.source == "AH") {
      if (verbose) message(
        "Loading ", REtype, " annotation data from RepeatMasker (hg38) database from AnnotationHub: ",
        remp_options(".default.AH.repeatmasker.hg38")
      )
      ah <- .initiateAH()
      if(is.null(ah)) {
        warning("AnnotationHub is currently not accessible, Please try again later.")
      } else {
        rmsk <- suppressMessages(ah[[remp_options(".default.AH.repeatmasker.hg38")]])
      }
    }
  }

  if (REtype == "Alu") {
    REFamily_grep <- remp_options(".default.AluFamily.grep")
  }
  if (REtype == "L1") {
    REFamily_grep <- remp_options(".default.L1Family.grep")
  }
  if (REtype == "ERV") {
    REFamily_grep <- remp_options(".default.ERVFamily.grep")
  }
  
  RE <- rmsk[grep(REFamily_grep, rmsk$repFamily)]
  RE <- RE[as.character(seqnames(RE)) %in% remp_options(".default.chr")] # chr1 - 22, chrX, chrY
  seqlevels(RE) <- remp_options(".default.chr") # remove uncommon chr

  RE <- sort(RE) # sort so that AH and UCSC are identical
  
  mcols(RE)$Index <- Rle(paste(REtype, formatC(seq_len(length(RE)),
    width = 7,
    flag = "0"
  ), sep = "_")) # Add internal index, meaningless for external database
  
  if(verbose) message("Obtained ", length(RE), " ", REtype, " (", genome, ").")
  
  return(RE)
}




#' @title Get RefSeq gene database
#'
#' @description
#' \code{fetchRefSeqGene} is used to obtain refSeq gene database provided by AnnotationHub (hg19) 
#' or UCSC web database (hg19/hg38).
#'
#' @param annotation.source Character parameter. Specify the source of annotation databases, including
#' the RefSeq Gene annotation database and RepeatMasker annotation database. If \code{"AH"}, the database 
#' will be obtained from the AnnotationHub package. If \code{"UCSC"}, the database will be downloaded 
#' from the UCSC website http://hgdownload.cse.ucsc.edu/goldenpath. The corresponding build (\code{"hg19"} or 
#' \code{"hg38"}) will be specified in the parameter \code{genome}.
#' @param genome Character parameter. Specify the build of human genome. Can be either \code{"hg19"} or 
#' \code{"hg38"}. Note that if \code{annotation.source == "AH"}, only hg19 database is available.
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
#' both \code{\link{GRanges}} object (for main refgene data) and \code{\link{GRangesList}} object
#' (for gene regions data).
#'
#' @examples
#' if (!exists("refgene.hg19")) 
#' refgene.hg19 <- fetchRefSeqGene(annotation.source = "AH", 
#'                                 genome = "hg19", 
#'                                 verbose = TRUE)
#' refgene.hg19
#' 
#' @export
fetchRefSeqGene <- function(annotation.source = c("AH", "UCSC"), 
                            genome = c("hg19", "hg38"), 
                            mainOnly = FALSE, 
                            verbose = FALSE) {
  annotation.source = match.arg(annotation.source)
  genome = match.arg(genome)
  
  if(genome == "hg19") {
    if(annotation.source == "AH") {
      if(verbose) message(
        "Loading static refSeq gene (hg19) database from AnnotationHub: ",
        remp_options(".default.AH.refgene.hg19"))
      ah <- .initiateAH()
      if(is.null(ah)) {
        message("AnnotationHub is currently not accessible, Switching to 'UCSC' source instead...")
        annotation.source <- "UCSC"
      } else refgene <- suppressMessages(ah[[remp_options(".default.AH.refgene.hg19")]])
    }
    if(annotation.source == "UCSC") {
      if(verbose) message("Loading UCSC refSeq gene (hg19) database from ", remp_options(".default.refGene.hg19.URL"))
      refgene <- .UCSCrefGeneDownload(url = remp_options(".default.refGene.hg19.URL"), 
                                      tag = "refGene.hg19",
                                      verbose)
    }
  }
  
  if(genome == "hg38") {
    if(annotation.source == "AH") {
      message("RefSeq Gene database in hg38 build is not available in AnnotationHub. Switching to 'UCSC' source instead...")
      annotation.source <- "UCSC"
    }
    if(annotation.source == "UCSC") {
      if(verbose) message("Loading UCSC refSeq gene (hg38) database from ", remp_options(".default.refGene.hg38.URL"))
      refgene <- .UCSCrefGeneDownload(url = remp_options(".default.refGene.hg38.URL"), 
                                      tag = "refGene.hg38",
                                      verbose)
    }
  }
  
  refgene <- refgene[as.character(seqnames(refgene)) %in%
    paste0("chr", c(seq_len(22), "X", "Y"))] ## chr1 - 22, chrX, chrY
  refgene$type <- substring(refgene$name, 1, 2) ## Protein coding gene (NM) or noncoding RNA gene (NR)
  refgene$index <- as.character(seq_len(length(refgene))) ## This is internal index, meaningless for external database

  ## Master database with detailed gene information
  refgene_main <- GRanges(
    seqnames = seqnames(refgene), ranges = ranges(refgene),
    strand = strand(refgene), name = mcols(refgene)$name, type = mcols(refgene)$type,
    index = mcols(refgene)$index
  )

  refseq2EG <- unlist(BiocGenerics::mget(refgene_main$name, org.Hs.eg.db::org.Hs.egREFSEQ2EG,
    ifnotfound = NA
  ))
  if (verbose) {
    message("Note: There are ", sum(is.na(refseq2EG)), " refseq ID cannot be mapped to Entrez ID, which could be obsolete.")
  }
  # obsoleteID <- refgene_main$name[is.na(refseq2EG)]
  refseq2EG <- na.omit(refseq2EG)

  EG2Symbol <- unlist(BiocGenerics::mget(refseq2EG, org.Hs.eg.db::org.Hs.egSYMBOL)) ## If EG -> symbol not found, this is a fatel error, should stop here

  refseq2EG_Symbol.DF <- DataFrame(
    refseqID = names(refseq2EG), EntrezGeneID = refseq2EG,
    GeneSymbol = EG2Symbol
  )

  refgene_main$EntrezGene <- refseq2EG_Symbol.DF$EntrezGeneID[base::match(
    refgene_main$name,
    refseq2EG_Symbol.DF$refseqID
  )]
  refgene_main$GeneSymbol <- refseq2EG_Symbol.DF$GeneSymbol[base::match(
    refgene_main$name,
    refseq2EG_Symbol.DF$refseqID
  )]

  if (mainOnly) {
    return(refgene_main)
  }

  ## TSS2000
  refgene_tss <- promoters(refgene,
    upstream = remp_options(".default.TSS.upstream"),
    downstream = remp_options(".default.TSS.downstream")
  )[, "index"]

  ## CDS
  refgene_cds <- GRanges(
    seqnames = seqnames(refgene), ranges = refgene$thick,
    strand = strand(refgene), index = mcols(refgene)$index
  )

  ## EXON
  block_shift <- shift(refgene$blocks, start(refgene) - 1) # minus one
  refgene_exon <- makeGRangesListFromFeatureFragments(
    seqnames = seqnames(refgene_main),
    fragmentStarts = start(block_shift), fragmentEnds = end(block_shift),
    strand = strand(refgene_main)
  )
  names(refgene_exon) <- refgene_main$index
  refgene_exon <- unlist(refgene_exon)
  refgene_exon$index <- as.integer(names(refgene_exon))
  names(refgene_exon) <- NULL

  ## UTR (5'UTR/3'UTR concept only applied for protein-coding gene)
  refgene_NM <- refgene[refgene$type == "NM", c("thick", "index")]
  refgene_f <- refgene_NM[strand(refgene_NM) == "+"]
  refgene_r <- refgene_NM[strand(refgene_NM) == "-"]

  refgene_5UTR_f <- GRanges(
    seqnames = seqnames(refgene_f),
    ranges = IRanges(
      start = start(refgene_f),
      end = start(refgene_f$thick)
    ),
    strand = strand(refgene_f),
    index = refgene_f$index
  )

  refgene_5UTR_r <- GRanges(
    seqnames = seqnames(refgene_r),
    ranges = IRanges(
      start = end(refgene_r$thick),
      end = end(refgene_r)
    ),
    strand = strand(refgene_r),
    index = refgene_r$index
  )

  refgene_3UTR_f <- GRanges(
    seqnames = seqnames(refgene_f),
    ranges = IRanges(
      start = end(refgene_f$thick),
      end = end(refgene_f)
    ),
    strand = strand(refgene_f),
    index = refgene_f$index
  )

  refgene_3UTR_r <- GRanges(
    seqnames = seqnames(refgene_r),
    ranges = IRanges(
      start = start(refgene_r),
      end = start(refgene_r$thick)
    ),
    strand = strand(refgene_r),
    index = refgene_r$index
  )

  refgene_5UTR <- c(refgene_5UTR_f, refgene_5UTR_r)
  refgene_3UTR <- c(refgene_3UTR_f, refgene_3UTR_r)

  return(list(
    main = refgene_main,
    regions = GRangesList(
      tss = refgene_tss,
      cds = refgene_cds,
      exon = refgene_exon,
      fiveUTR = refgene_5UTR,
      threeUTR = refgene_3UTR
    )
  ))
}



#' @title Find RE-CpG genomic location given RE ranges information
#'
#' @description
#' \code{findRECpG} is used to obtain RE-CpG genomic location data.
#'
#' @param RE A \code{\link{GRanges}} object of RE genomic location database. This
#' can be obtained by \code{\link{fetchRMSK}}.
#' @param REtype Type of RE. Currently \code{"Alu"}, \code{"L1"}, and \code{"ERV"} are supported.
#' @param genome Character parameter. Specify the build of human genome. Can be either \code{"hg19"} or 
#' \code{"hg38"}. User should make sure the genome build of \code{RE} is consistent with this parameter.
#' @param be A \code{\link{BiocParallel}} object containing back-end information that is
#' ready for parallel computing. This can be obtained by \code{\link{getBackend}}. If not specified,
#' non-parallel mode is used.
#' @param verbose logical parameter. Should the function be verbose?
#'
#' @details
#' CpG site is defined as 5'-C-p-G-3'. It is reasonable to assume that the methylation
#' status across all CpG/CpG dyads are concordant. Maintenance methyltransferase exhibits
#' a preference for hemimethylated CpG/CpG dyads (methylated on one strand only).
#' As a result, methyaltion status of CpG sites in both forward and reverse strands are usually consistent.
#' Therefore, to accommodate the cytosine loci in both strands, the returned genomic
#' ranges cover the 'CG' sequence with width of 2. The 'strand' information indicates the strand of the RE.
#' Locating CpG sites in RE sequences can be computation intensive. It is recommanded to get more than
#' one work in the backend for a faster running speed.
#'
#' @return A \code{\link{GRanges}} object containing identified RE-CpG genomic
#' location data.
#'
#' @examples
#' data(Alu.hg19.demo)
#' RE.CpG <- findRECpG(RE = Alu.hg19.demo, 
#'                     REtype = "Alu", 
#'                     genome = "hg19", 
#'                     verbose = TRUE)
#' RE.CpG
#' 
#' @export
findRECpG <- function(RE, 
                      REtype = c("Alu", "L1", "ERV"), 
                      genome = c("hg19", "hg38"), 
                      be = NULL, 
                      verbose = FALSE) {
  if (is.null(be)) be <- getBackend(1, verbose)
  REtype <- match.arg(REtype)
  genome <- match.arg(genome)
  
  ## Flank RE by 1 base pair to cover the cases where CG is located in
  ## RE's head or tail
  RE <- .twoWayFlank(RE, 1)

  if (verbose) {
    message(
      "    Getting sequence of the total ", length(RE),
      " ", REtype, " ..."
    )
  }
  
  if(genome == "hg19") {
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
      if(isNamespaceLoaded("BSgenome.Hsapiens.UCSC.hg38")) unloadNamespace("BSgenome.Hsapiens.UCSC.hg38")
      if(!isNamespaceLoaded("BSgenome.Hsapiens.UCSC.hg19")) attachNamespace("BSgenome.Hsapiens.UCSC.hg19")
      if(verbose) message("    Using UCSC Human Reference Genome: hg19")
      hs <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    } else stop("Please install missing package: BSgenome.Hsapiens.UCSC.hg19") 
  }
  
  if(genome == "hg38") {
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      if(isNamespaceLoaded("BSgenome.Hsapiens.UCSC.hg19")) unloadNamespace("BSgenome.Hsapiens.UCSC.hg19")
      if(!isNamespaceLoaded("BSgenome.Hsapiens.UCSC.hg38")) attachNamespace("BSgenome.Hsapiens.UCSC.hg38")
      if(verbose) message("    Using UCSC Human Reference Genome: hg38")
      hs <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    } else stop("Please install missing package: BSgenome.Hsapiens.UCSC.hg38") 
  }
  
  SEQ.RE <- BSgenome::getSeq(hs, 
                             GRanges(seqnames(RE), 
                                     IRanges(start = start(RE), end = end(RE)))) # ignore strand

  if (verbose) {
    message("    Identifying CpG sites in ", REtype, " sequence ...")
  }

  if (BiocParallel::bpnworkers(be) > 1) {
    BiocParallel::bpstart(be)
    .bploadLibraryQuiet("Biostrings", be)
    RE.CpG <- bpvec(SEQ.RE, .vRECpGPos, CpG = DNAString("CG"), BPPARAM = be)
    BiocParallel::bpstop(be)
  } else {
    RE.CpG <- .vRECpGPos(SEQ.RE, CpG = DNAString("CG"))
  }

  noCpGind <- which(vapply(RE.CpG, length, integer(1)) == 0) ## Remove sequence without CpG sites

  if (length(noCpGind) == 0) {
    if (verbose) {
      message("    All ", REtype, " sequences contain CpGs.")
    }
  } else {
    if (verbose) {
      message("    There are ", length(noCpGind), " ", REtype, " sequences contain no CpG site.")
    }
    RE <- RE[-noCpGind]
    RE.CpG <- RE.CpG[-noCpGind]
  }

  ## Extend 'C' to 'CG' AND Shift to the real genomic location
  RE.CpG.IRList <- IRangesList(start = RE.CpG - 1, end = RE.CpG)
  block_shift <- shift(RE.CpG.IRList, start(RE))

  RE.CpG <- makeGRangesListFromFeatureFragments(
    seqnames = seqnames(RE),
    fragmentStarts = start(block_shift),
    fragmentEnds = end(block_shift),
    strand = strand(RE)
  )
  elementLen <- elementNROWS(RE.CpG)

  RE.CpG <- unlist(RE.CpG)
  RE.CpG$Index <- Rle(runValue(RE$Index), elementLen)
  RE.CpG <- RE.CpG[order(RE.CpG$Index)]

  if (verbose) {
    message(
      "    Identified ", length(RE.CpG), " CpG/CpG dyads in ",
      nrun(RE.CpG$Index), " ", REtype
    )
  }

  return(RE.CpG)
}



#' @title Annotate genomic ranges data with gene region information.
#'
#' @description
#' \code{GRannot} is used to annotate a \code{\link{GRanges}} dataset with gene region
#' information using refseq gene database
#'
#' @param object.GR An \code{\link{GRanges}} object of a genomic location database.
#' @param refgene A complete refGene annotation database returned by
#' \code{\link{fetchRefSeqGene}} (with parameter \code{mainOnly = FALSE}).
#' @param symbol Logical parameter. Should the annotation return gene symbol?
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
#' data(Alu.hg19.demo)
#' if (!exists("refgene.hg19")) 
#'   refgene.hg19 <- fetchRefSeqGene(annotation.source = "AH", 
#'                                   genome = "hg19",
#'                                   verbose = TRUE)
#' Alu.hg19.demo.refGene <- GRannot(Alu.hg19.demo, refgene.hg19, verbose = TRUE)
#' Alu.hg19.demo.refGene
#' 
#' @export
GRannot <- function(object.GR, 
                    refgene, 
                    symbol = FALSE, 
                    verbose = FALSE) {

  ## Check if the object.GR contains a metadata column called "Index"
  if (!"Index" %in% colnames(mcols(object.GR))) {
    object.GR$Index <- as.character(seq_len(length(object.GR)))
  }

  if (any(duplicated(as.character(object.GR$Index)))) {
    stop("The 'Index' column provided in the GRanges object must be unique.")
  }

  main <- refgene$main

  ####### InNM
  NM.GR <- main[main$type == "NM"]
  NM.GR <- .oneWayFlank(NM.GR, 2000, start = TRUE)
  object.GR <- .addannot(
    object.GR,
    NM.GR,
    main,
    symbol,
    "InNM",
    verbose
  )

  ####### InNR
  NR.GR <- main[main$type == "NR"]
  NR.GR <- .oneWayFlank(NR.GR, 2000, start = TRUE)
  object.GR <- .addannot(
    object.GR,
    NR.GR,
    main,
    symbol,
    "InNR",
    verbose
  )

  ####### InTSS2000
  object.GR <- .addannot(
    object.GR,
    refgene$regions$tss,
    main,
    symbol,
    "InTSS",
    verbose
  )

  ####### In5UTR
  object.GR <- .addannot(
    object.GR,
    refgene$regions$fiveUTR,
    main,
    symbol,
    "In5UTR",
    verbose
  )

  ####### InCDS
  object.GR <- .addannot(
    object.GR,
    refgene$regions$cds,
    main,
    symbol,
    "InCDS",
    verbose
  )

  ####### InExon
  object.GR <- .addannot(
    object.GR,
    refgene$regions$exon,
    main,
    symbol,
    "InExon",
    verbose
  )

  ####### In3UTR
  object.GR <- .addannot(
    object.GR,
    refgene$regions$threeUTR,
    main,
    symbol,
    "In3UTR",
    verbose
  )

  return(object.GR)
}



## --------------------------------------
## Internal functions

# Used by findRECpG (vetorized function)
.vRECpGPos <- function(seq, CpG) {
  start(vmatchPattern(CpG, seq))
}

# Used by addannot
.clps <- function(x) paste0(x, collapse = "|")

# Used by GRannot
.addannot <- function(object.GR, annot.GR, main, symbol, annotName, verbose) {
  if (verbose) {
    message("    ", annotName, " ... ", appendLF = FALSE)
  }

  object_annot_Hit <- findOverlaps(object.GR, annot.GR, ignore.strand = TRUE)
  object_raw <- DataFrame(InRegion = annot.GR[subjectHits(object_annot_Hit)]$index)
  InRegion <- rep(NA, length(object.GR))

  if (nrow(object_raw) > 0) {
    if (symbol) {
      object_raw$InRegion <- main$GeneSymbol[as.integer(object_raw$InRegion)]
    }

    object_agg <- aggregate(
      object_raw,
      list(Index = as.character(object.GR$Index[queryHits(object_annot_Hit)])),
      .clps
    )
    InRegion[base::match(object_agg$Index, 
                         runValue(object.GR$Index))] <- object_agg$InRegion
  }

  object.GR$newRegion <- InRegion
  mcols(object.GR) <- .changeColNames(
    mcols(object.GR),
    "newRegion",
    annotName
  )

  if (verbose) {
    message("  ", nrow(object_agg), " found!")
  }

  return(object.GR)
}

# Used by fetchRefSeqGene
.UCSCrefGeneDownload <- function(url, tag, verbose)
{
  refgene_raw <- .webDownload(url = url, 
                              tag = tag, 
                              col_types = readr::cols_only(X2 = readr::col_character(),
                                                           X3 = readr::col_character(),
                                                           X4 = readr::col_character(),
                                                           X5 = readr::col_integer(),
                                                           X6 = readr::col_integer(),
                                                           X7 = readr::col_integer(),
                                                           X8 = readr::col_integer(),
                                                           X10 = readr::col_character(),
                                                           X11 = readr::col_character()), 
                              verbose = verbose)
  refgene <- GRanges(seqnames = refgene_raw$X3, IRanges(start = refgene_raw$X5 + 1, end = refgene_raw$X6),
                     strand = refgene_raw$X4, name = vapply(strsplit(refgene_raw$X2, "\\."), function(x) x[1], character(1)), 
                     thick = IRanges(start = refgene_raw$X7 + 1, end = refgene_raw$X8),
                     blocks = IRangesList(start = mapply(function(x, y) x + 1 - y,  
                                                         lapply(strsplit(refgene_raw$X10, ","), as.integer), 
                                                         refgene_raw$X5), 
                                          end = mapply(function(x, y) x - y,  
                                                       lapply(strsplit(refgene_raw$X11, ","), as.integer), 
                                                       refgene_raw$X5)))
  return(refgene)
}
