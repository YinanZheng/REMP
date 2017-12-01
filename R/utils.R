## Internal functions:

.timeTrace <- function(startTime, indent = "  ") {
  t <- round(as.numeric(Sys.time() - startTime, units = "secs"), 0)
  paste0(indent, "(", t, " sec.)")
}

.forwardSlashPath <- function(path) {
  gsub("\\\\", "/", path)
}

# Assign the pacakge to workers and suppress the loading messsage!
.bploadLibraryQuiet <- function(libraryList, be) {
  noreturn <- BiocParallel::bplapply(seq_len(BiocParallel::bpnworkers(be)), function(i, libraryList) {
    suppressPackageStartupMessages(lapply(libraryList, library, character.only = TRUE))
  }, BPPARAM = be, libraryList = libraryList)
}

.iblkrow <- function(M, chunks) {
  i <- 1
  uni_row <- nrow(M)
  it <- idiv(uni_row, chunks = chunks)
  
  nextEl <- function() {
    if (i == uni_row + 1) 
      return(NULL)
    n <- nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[r, , drop = FALSE]
  }
  nextEl
}

.iblkrow_dup <- function(M, chunks, index_col_name) {
  i <- 1
  ind <- which(colnames(M) == index_col_name)
  index <- M[,ind]
  unique_index <- unique(index)
  unique_length <- length(unique_index)
  
  it <- idiv(unique_length, chunks = chunks)
  
  nextEl <- function() {
    if (i == unique_length + 1) 
      return(NULL)
    n <- nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[index %in% unique_index[r], , drop = FALSE]
  }
  nextEl
}

.aggregateREMP <- function(dat, index_col_name)
{
  eval(parse(text = paste0("aggregate(.~",index_col_name,", dat, mean, na.rm = TRUE, na.action = na.pass)")))
}

.guessArrayType <- function(methyDat) {
  nProbes <- nrow(methyDat)
  if (nProbes <= remp_options(".default.27k.total.probes"))
    return("27k")
  if (nProbes > remp_options(".default.27k.total.probes") & 
      nProbes <= remp_options(".default.450k.total.probes"))
    return("450k")
  if (nProbes > remp_options(".default.450k.total.probes") & 
      nProbes <= remp_options(".default.epic.total.probes")) 
    return("EPIC")
  if (nProbes > remp_options(".default.epic.total.probes")) 
    return("UNKNOWN")
}

.guessBetaorM <- function(methyDat) {
  methyDat_sample <- methyDat[sample(seq_len(nrow(methyDat)), min(5000, nrow(methyDat))), 
                              sample(seq_len(ncol(methyDat)), min(5, ncol(methyDat))), drop = FALSE]
  rng <- range(methyDat_sample, na.rm = TRUE)
  if (rng[1] >= 0 & rng[2] <= 1 | 
      sum(methyDat_sample < 0 | methyDat_sample > 1, na.rm = TRUE) < 
      nrow(methyDat_sample) * ncol(methyDat_sample) * 0.01) 
    return("beta") else return("M")
}

.twoWayFlank <- function(object.GR, width) {
  .isGROrStop(object.GR)
  start(object.GR) <- start(object.GR) - width
  end(object.GR) <- end(object.GR) + width
  return(object.GR)
}

.oneWayFlank <- function(object.GR, width, start = TRUE) {
  .isGROrStop(object.GR)
  if(any(runValue(strand(object.GR)) == "*"))
    stop("Strand information must not be missing.")
  forward_ind <- strand(object.GR)=="+"
  reverse_ind <- strand(object.GR)=="-"
  
  if(start)
  {
    start(object.GR[forward_ind]) <- start(object.GR[forward_ind]) - width
    end(object.GR[reverse_ind]) <- end(object.GR[reverse_ind]) + width
  } else {
    end(object.GR[forward_ind]) <- end(object.GR[forward_ind]) + width
    start(object.GR[reverse_ind]) <- start(object.GR[reverse_ind]) - width
  }
  return(object.GR)
}

.changeColNames <- function(DForGR, oldnames, newnames) {
  if(is(DForGR, "DataFrame"))
  {
    cnames <- colnames(DForGR)
    colnames(DForGR)[cnames %in% oldnames] <- newnames
    return(DForGR)
  }
  if(is(DForGR, "GRanges"))
  {
    cnames <- colnames(mcols(DForGR))
    colnames(mcols(DForGR))[cnames %in% oldnames] <- newnames
    return(DForGR)
  }
}

#### Message display function

.showREMParceInfo <- function(object)
{
  info <- object@REMParcelInfo
  cat("REMParcel object\n")
  cat("RE type:", info[["REtype"]], "\n")
  cat("Illumina platform:", info[["platform"]], "\n")
  cat("Valid (max)", paste0(info[["REtype"]], "-CpG"), "flanking window size:", info[["max.win"]], "\n")
  cat("Number of RE:", length(object@RE), "\n")
  cat("Number of", paste0(info[["REtype"]], "-CpG:"), length(object@RECpG), "\n")
}

.showTrainingStats <- function(REStats, REtype, printType, indent) {
  if(!isEmpty(REStats))
  {
    p <- rep(NA, 2)
    
    p[1] <- paste0(indent, REStats[1, 1], " profiled ", 
                   REtype, " by Illumina array are used for model training.")
    p[2] <- paste0(indent, REStats[1, 2], 
                   " ", REtype, "-CpGs that have at least 2 neighboring profiled CpGs are used for model training.")
  } else {
    p <- "N/A"
  } 
  
  if(printType == "message")
    trashbin <- sapply(p, message)
  if(printType == "cat")
    trashbin <- sapply(p, function(x) cat(x, "\n"))
}

.showREStats <- function(REStats, REtype, printType, indent, notAggregated) {
  if(!isEmpty(REStats))
  {
    p <- paste0(indent, "The data cover ", REStats[1, 3], " ", REtype)
    if(notAggregated) p <- paste0(p, " (", REStats[1, 4], " ", REtype, "-CpG).")
  } else {
    p <- "N/A"
  } 
  
  if(printType == "message")
    trashbin <- sapply(p, message)
  if(printType == "cat")
    trashbin <- sapply(p, function(x) cat(x, "\n"))
}

.showGeneStats <- function(GeneStats, REtype, printType, indent) {
  if(!isEmpty(GeneStats))
  {
    p <- rep(NA, 4)
    
    p[1] <- paste0(indent, "Gene coverage by ", REtype, " (out of total refSeq Gene):")
    p[2] <- paste0(indent, indent, GeneStats[2, 3], 
                   " (", round(GeneStats[2, 3]/GeneStats[1, 3] * 100, 2), "%) total genes;")
    p[3] <- paste0(indent, indent, GeneStats[2, 1], 
                   " (", round(GeneStats[2, 1]/GeneStats[1, 1] * 100, 2), "%) protein-coding genes;")
    p[4] <- paste0(indent, indent, GeneStats[2, 2], 
                   " (", round(GeneStats[2, 2]/GeneStats[1, 2] * 100, 2), "%) non-coding RNA genes.")
  } else {
    p <- "N/A"
  } 
  
  if(printType == "message")
    trashbin <- sapply(p, message)
  if(printType == "cat")
    trashbin <- sapply(p, function(x) cat(x, "\n"))
}

.showSampleID <- function(object, indent) {
  x <- object@sampleID
  cat("Number of sample(s):", length(x), "\n")
  cat("Sample ID:")
  if (length(x) > 5) {
    x <- x[seq_len(5)]
    cat(indent, paste0(x, collapse = ", "), " ...\n")
  } else {
    cat(indent, paste0(x, collapse = ", "), "\n")
  }
}

.showREMPinfo <- function(object) {
  info <- object@REMPInfo
  cat("RE type:", info[["REtype"]], "\n")
  cat("Methylation profiling platform:", info[["platform"]], "\n")
  cat("Flanking window size:", info[["win"]], "\n")
  cat("Prediction model:", info[["predictModel"]], "\n")
  cat("QC model:", info[["QCModel"]], "\n")
  if(grepl("Random Forest", info[["predictModel"]]))
  {
    cat("Seed:", metadata(object)$Seed, "\n")
  }  
  restats <- metadata(object)$REStats
  if(!isEmpty(restats))
  {
    if(grepl("aggregated", info[["REtype"]]))
    {
      cat("Covered", restats[1,3], info[["REtype"]], "\n")
    } else {
      cat("Covered", restats[1,4], "CpG sites in", 
          restats[1,3], info[["REtype"]], "\n")
    }
  }
}

.showCpGcountbyChr <- function(object){
  info <- object@REMPInfo
  if(grepl("aggregated", info[["REtype"]]))
  {
    cat("Number of", info[["REtype"]], "by chromosome:")
  } else {
    cat("Number of", paste0(info[["REtype"]],"-CpGs"), "by chromosome:")
  }
  chr <- table(as.character(seqnames(rowRanges(object))))
  chr.name <- names(chr)
  chr.name[chr.name=="chrX"] = "chr23"
  chr.name[chr.name=="chrY"] = "chr24"
  chr.num <- as.numeric(substring(chr.name, 4, 5))
  chr <- chr[order(chr.num)]
  
  for(i in seq_len(ceiling(length(chr) / 8)))
  {
    print(chr[((i-1) * 8 + 1) : min((8 * i), length(chr))])
  }
}

.showPredictionSummary <- function(object){
  options(digits=10)
  Bvalue <- assays(object)[["rempB"]]
  cat("Distribution of methylation value (beta value):", "\n")
  print(summary(as.numeric(Bvalue)))
  options(digits=7)
}

.showQCSummary <- function(object){
  options(digits=10)
  QC <- assays(object)[["rempQC"]]
  cat("Distribution of reliability score:", "\n")
  print(summary(as.numeric(QC)))
  options(digits=7)
}



## Convert refgene index to gene symbol or entrez id

.decodeEntrez <- function(x, refgene) {
  paste0(unique(refgene[x, "EntrezGene"]), collapse = "|")
}

.decodeSymbol <- function(x, refgene) {
  paste0(unique(refgene[x, "GeneSymbol"]), collapse = "|")
}



## Methylation beta <--> M value conversion

.toBeta <- function(M) {
  2^(M)/(1 + 2^(M))
}

.toM <- function(beta) {
  log2(beta / (1 - beta))
}



## Class check

.isGROrStop <- function(object) {
  if (!is(object, "GRanges")) 
    stop(sprintf("object is of class '%s', but needs to be of class 'GRanges' (see object definition in package 'GenomicRanges')", 
                 class(object)))
}

.isRSetOrStop <- function(object) {
  if (!is(object, "RatioSet") && !is(object, "GenomicRatioSet")) 
    stop(sprintf("object is of class '%s', but needs to be of class '[Genomic]RatioSet' (see object definition in package 'minfi')", 
                 class(object)))
}

.isREMParcelOrStop <- function(object) {
  if (!is(object, "REMParcel")) 
    stop(sprintf("object is of class '%s', but needs to be of class 'REMParcel'.", 
                 class(object)))
}

.isREMProductOrStop <- function(object) {
  if (!is(object, "REMProduct")) 
    stop(sprintf("object is of class '%s', but needs to be of class 'REMProduct'.", 
                 class(object)))
}

