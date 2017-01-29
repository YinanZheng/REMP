
## REMProduct methods

#### Accessors

#' @rdname REMProduct-class
setMethod("rempM", signature(object = "REMProduct"), function(object) {
  M <- assays(object)[["rempM"]]
  return(DataFrame(M, check.names = FALSE))
})
# rempM(object)

#' @rdname REMProduct-class
setMethod("rempB", signature(object = "REMProduct"), function(object) {
  B <- assays(object)[["rempB"]]
  return(DataFrame(B, check.names = FALSE))
})
# rempB(object)

#' @rdname REMProduct-class
setMethod("rempQC", signature(object = "REMProduct"), function(object) {
  QC <- assays(object)[["rempQC"]]
  return(DataFrame(QC, check.names = FALSE))
})
# rempQC(object)

#' @rdname REMProduct-class
setMethod("imp", signature(object = "REMProduct"), function(object) {
  return(metadata(object)$varImp)
})
# imp(object)

#' @rdname REMProduct-class
setMethod("annotation", signature(object = "REMProduct"), function(object) {
  return(metadata(object)$REannotation)
})
# annotation(object)

#' @rdname REMProduct-class
setMethod("stats", signature(object = "REMProduct"), function(object) {
  return(list(RE_Statistics = metadata(object)$REStats,
              Gene_Statistics = metadata(object)$GeneStats))
})
# stats(object)


#### Utilities
#' @rdname REMProduct-class
setMethod("plot", signature(x = "REMProduct", y = "missing"),
          function(x, type = c("individual", "overall"), ...) {
  type <- match.arg(type)
  default_main <- ""
  default_xlab <- "Methylation value (beta)"
  dotdotdot <- list(...)
  if(!hasArg(main))
    dotdotdot$main <- default_main
  if(!hasArg(xlab))
    dotdotdot$xlab <- default_xlab

  Bval <- assays(x)[["rempB"]]

  if(type == "individual")
  {
    do.call(plot, c(list(x = density(Bval[,1])), dotdotdot))
    if(ncol(Bval) >= 2){
      for(i in seq(2, ncol(Bval)))
      {
        do.call(lines, c(list(x = density(Bval[,i])), dotdotdot))
      }
    }
  }

  if(type == "overall")
  {
    Bval.mean = rowMeans(Bval, na.rm = TRUE)
    do.call(plot, c(list(x = density(Bval.mean)), dotdotdot))
  }
})
# plot(object)
# plot(object, type = "overall")

#' @rdname REMProduct-class
setMethod("details", signature(object = "REMProduct"), function(object) {
  .showREMPinfo(object)
  cat("\n")
  
  .showCpGcountbyChr(object)
  cat("\n")
  
  cat("Coverage information:\n")
  .showREStats(metadata(object)$REStats, object@REMPInfo[["REtype"]], "cat", "  ")
  .showGeneStats(metadata(object)$GeneStats, object@REMPInfo[["REtype"]], "cat", "  ")
  cat("\n")
  
  .showPredictionSummary(object)
  if (object@REMPInfo[["QCModel"]] != "N/A") {
    cat("\n")
    .showQCSummary(object)
  }
})
# details(object)

#' @rdname REMProduct-class
setMethod("decodeAnnot", signature(object = "REMProduct"), 
          function(object, type = c("symbol", "entrez")) {
  .isREMProductOrStop(object)
  skip <- TRUE  
  type <- match.arg(type)
  message("Decoding ", object@REMPInfo[["REtype"]], " annotation to ", type, " ...")
  
  annot <- metadata(object)$REannotation
  code <-  metadata(object)$regionCode
  refgene <- metadata(object)$refGene
  containSymbol <- any(grepl("symbol", colnames(mcols(annot))))
  containEntrez <- any(grepl("entrez", colnames(mcols(annot))))
    
  if(type == "symbol" & !containSymbol)
  {
    skip = FALSE
    decodeFun <- .decodeSymbol
  }
  if(type == "entrez" & !containEntrez)
  {
    skip = FALSE
    decodeFun <- .decodeEntrez
  }
  
  if(!skip)
  {
    for (region in remp_options(".default.genomicRegionColNames")) {
      InRegion.original <- code[, region]
      ind.missing <- which(is.na(InRegion.original))
      InRegion <- na.omit(InRegion.original)
      res <- strsplit(InRegion, "[|]")
      res <- lapply(res, as.numeric)
      res <- sapply(res, decodeFun, refgene = refgene)
      InRegion.original[-ind.missing] <- res
      annot$newregionEncode <- InRegion.original
      annot <- .changeColNames(annot, "newregionEncode", 
                               paste(region, type, sep = "."))
    }
    metadata(object)$REannotation <- annot
  }
  return(object)
})
# annotation(object)
# remptest <- decodeAnnot(object)
# annotation(remptest)
# remptest <- decodeAnnot(remptest)
# annotation(remptest)
# remptest <- decodeAnnot(remptest, type = "entrez")
# annotation(remptest)

#' @rdname REMProduct-class
setMethod("trim", signature(object = "REMProduct"), 
          function(object, threshold = 1.7, missingRate = 0.2) {
            method <- object@REMPInfo[["predictModel"]]
            if(method != "Random Forest")
            {
              message("Trim is only applicable to prediction using Random Forest model. No changes made.")
              return(object)
            }
            
            REtype = object@REMPInfo[["REtype"]]
            beta <- as.matrix(rempB(object))
            M <- as.matrix(rempB(object))
            QC <- as.matrix(rempQC(object))
            badInd <- QC > threshold
            QC[badInd] <- NA
            beta[badInd] <- NA
            M[badInd] <- NA
            
            removeCpG <- rowMeans(badInd) > missingRate
            
            RE_annotation <- annotation(object)
            RE_CpG_ILMN <-  metadata(object)$RECpG
            regionCode <-  metadata(object)$regionCode
            refgene_main <- metadata(object)$refGene

            beta <- beta[!removeCpG,,drop = FALSE]
            M <- M[!removeCpG,,drop = FALSE]
            QC <- QC[!removeCpG,,drop = FALSE]
            cpgRanges <- rowRanges(object)[!removeCpG, ]
            
            trimmed_RE_list <- runValue(cpgRanges$RE.Index)
            regionCode <- regionCode[RE_annotation$Index %in% trimmed_RE_list,]
            RE_annotation <- RE_annotation[RE_annotation$Index %in% trimmed_RE_list]
            
            ## Updated RE coverage
            RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN, 
                                             REtype, indent = "    ", TRUE)
            
            # Updated Gene coverage
            GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main, 
                                                 REtype, indent = "    ", TRUE)
            
            
            ## Update object
            object_trim <- REMProduct(REtype = object@REMPInfo[["REtype"]], 
                                      platform = object@REMPInfo[["platform"]], 
                                      win = object@REMPInfo[["win"]],
                                      predictModel = paste0(object@REMPInfo[["predictModel"]], " - trimmed (", threshold, ")"),  
                                      QCModel = object@REMPInfo[["QCModel"]], 
                                      rempM = M, rempB = beta, rempQC = QC,
                                      cpgRanges = cpgRanges, sampleInfo = colData(remp.res),
                                      REannotation = RE_annotation, 
                                      RECpG = RE_CpG_ILMN,
                                      regionCode = regionCode,
                                      refGene = refgene_main,
                                      varImp = imp(object), 
                                      REStats = RE_COVERAGE, GeneStats = GENE_COVERAGE)
            
            return(object_trim)
          })
# object_trim <- trim(object)
# details(object_trim)
