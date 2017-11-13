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
setMethod("rempImp", signature(object = "REMProduct"), function(object) {
  return(metadata(object)$varImp)
})
# rempImp(object)

#' @rdname REMProduct-class
setMethod("rempAnnot", signature(object = "REMProduct"), function(object) {
  return(metadata(object)$REannotation)
})
# rempAnnot(object)

#' @rdname REMProduct-class
setMethod("rempStats", signature(object = "REMProduct"), function(object) {
  return(list(RE_Statistics = metadata(object)$REStats,
              Gene_Statistics = metadata(object)$GeneStats))
})
# rempStats(object)


#### Utilities
#' @rdname REMProduct-class
setMethod("plot", signature(x = "REMProduct", y = "missing"),
          function(x, type = c("individual", "overall"), ...) {
  type <- match.arg(type)
  default_main <- ""
  default_xlab <- "Methylation value (beta)"
  default_xlim <- c(0,1)
  dotdotdot <- list(...)
  if(!hasArg("main"))
    dotdotdot$main <- default_main
  if(!hasArg("xlab"))
    dotdotdot$xlab <- default_xlab
  if(!hasArg("xlim"))
    dotdotdot$xlim <- default_xlim
    
  Bval <- rempB(x)

  if(type == "individual")
  {
    do.call(plot, c(list(x = density(na.omit(Bval[,1]))), dotdotdot))
    if(ncol(Bval) >= 2){
      for(i in seq(2, ncol(Bval)))
      {
        do.call(lines, c(list(x = density(na.omit(Bval[,i]))), dotdotdot))
      }
    }
  }

  if(type == "overall")
  {
    Bval.mean = rowMeans(Bval, na.rm = TRUE)
    do.call(plot, c(list(x = density(na.omit(Bval.mean))), dotdotdot))
  }
})
# plot(object)
# plot(object, type = "overall")

#' @rdname REMProduct-class
setMethod("details", signature(object = "REMProduct"), function(object) {
  REtype = object@REMPInfo[["REtype"]]
  notAggregated <- TRUE
  if(grepl("aggregated", REtype)) notAggregated <- FALSE
  
  method = object@REMPInfo[["predictModel"]]
  notTrimmed <- TRUE
  notNaive <- TRUE
  if(grepl("trimmed", method)) notTrimmed <- FALSE
  if(grepl("Naive", method)) notNaive <- FALSE
  
  .showREMPinfo(object)
  cat("\n")
  
  .showCpGcountbyChr(object)
  cat("\n")
  
  if(notAggregated & notTrimmed & notNaive)
  {
    cat("Training information:\n")
    .showTrainingStats(metadata(object)$REStats, object@REMPInfo[["REtype"]], "cat", "  ")
    cat("\n")
  }
  
  cat("Coverage information:\n")
  .showREStats(metadata(object)$REStats, object@REMPInfo[["REtype"]], "cat", "  ", notAggregated)
  .showGeneStats(metadata(object)$GeneStats, object@REMPInfo[["REtype"]], "cat", "  ")
  cat("\n")
  
  .showPredictionSummary(object)
  if (object@REMPInfo[["QCModel"]] != "N/A") {
    cat("\n")
    .showQCSummary(object)
  }
})

#' @rdname REMProduct-class
setMethod("decodeAnnot", signature(object = "REMProduct"), 
          function(object, type = c("symbol", "entrez"), ncore = NULL, BPPARAM = NULL) {
  .isREMProductOrStop(object)
            
  if (is.null(ncore)) 
    ncore <- parallel::detectCores()
  
  be <- getBackend(ncore, BPPARAM, TRUE)
  
  skip <- TRUE  
  type <- match.arg(type)
  message("Decoding ", object@REMPInfo[["REtype"]], " annotation to ", type, " ...")
  
  annot <- metadata(object)$REannotation
  code <-  metadata(object)$regionCode
  refgene <- as.data.frame(metadata(object)$refGene)
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
    bpstart(be)
    .bploadLibraryQuiet("IRanges", be)
    # message("Pacakge loaded!")
    for (region in remp_options(".default.genomicRegionColNames")) {
      InRegion.original <- code[, region]
      ind.missing <- which(is.na(InRegion.original))
      InRegion <- na.omit(InRegion.original)
      res <- strsplit(InRegion, "[|]")
      res <- lapply(res, as.numeric)
      res <- IRanges::NumericList(res)
      res <- unlist(bplapply(res, decodeFun, refgene = refgene, BPPARAM=be))
      InRegion.original[-ind.missing] <- res
      annot$newregionEncode <- InRegion.original
      annot <- .changeColNames(annot, "newregionEncode", 
                               paste(region, type, sep = "."))
    }
    metadata(object)$REannotation <- annot
    bpstop(be)
  }
  return(object)
})

#' @rdname REMProduct-class
setMethod("rempTrim", signature(object = "REMProduct"), 
          function(object, threshold = 1.7, missingRate = 0.2) {
            method <- object@REMPInfo[["predictModel"]]
            
            if(grepl("trimmed", method))
            {
              previous_thres <- gsub("[\\(\\)]", "", regmatches(method, gregexpr("\\(.*?\\)", method))[[1]])
              if (threshold >= as.numeric(previous_thres))
              {
                message("More stringent or equal threshold (", previous_thres, ") has been applied. No changes made.")
                return(object)
              } else {
                object@REMPInfo[["predictModel"]] = method <- "Random Forest"
              }
            }
            
            if(!grepl("Random Forest", method))
            {
              message("rempTrim() is only applicable to prediction using Random Forest model. No changes made.")
              return(object)
            }
            
            REtype = object@REMPInfo[["REtype"]]
            # beta <- as.matrix(rempB(object))
            M <- as.matrix(rempM(object))
            QC <- as.matrix(rempQC(object))
            badInd <- QC > threshold
            QC[badInd] <- NA
            # beta[badInd] <- NA
            M[badInd] <- NA
            
            removeCpG <- rowMeans(badInd) > missingRate
            
            RE_annotation <- rempAnnot(object)
            RE_CpG_ILMN <-  metadata(object)$RECpG
            regionCode <-  metadata(object)$regionCode
            refgene_main <- metadata(object)$refGene

            # beta <- beta[!removeCpG,,drop = FALSE]
            M <- M[!removeCpG,,drop = FALSE]
            QC <- QC[!removeCpG,,drop = FALSE]
            cpgRanges <- rowRanges(object)[!removeCpG, ]
            
            mcols(RE_annotation) <- cbind(mcols(RE_annotation), regionCode)
            RE_annotation <- subsetByOverlaps(RE_annotation, cpgRanges)
            RE_annotation_name <- colnames(mcols(RE_annotation))
            regionCode <- mcols(RE_annotation)[remp_options(".default.genomicRegionColNames")]
            RE_annotation <- RE_annotation[, RE_annotation_name[!RE_annotation_name %in% 
                                                                  remp_options(".default.genomicRegionColNames")]]
            
            ## Updated RE coverage
            RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN, 
                                             REtype, indent = "    ", FALSE)
            
            # Updated Gene coverage
            GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main, 
                                                 REtype, indent = "    ", FALSE)
            
            
            ## Update object
            object_trim <- REMProduct(REtype = object@REMPInfo[["REtype"]], 
                                      platform = object@REMPInfo[["platform"]], 
                                      win = object@REMPInfo[["win"]],
                                      predictModel = paste0(object@REMPInfo[["predictModel"]], " - trimmed (", threshold, ")"),  
                                      QCModel = object@REMPInfo[["QCModel"]], 
                                      rempM = M, rempQC = QC,
                                      cpgRanges = cpgRanges, sampleInfo = colData(object),
                                      REannotation = RE_annotation, 
                                      RECpG = RE_CpG_ILMN,
                                      regionCode = regionCode,
                                      refGene = refgene_main,
                                      varImp = rempImp(object), 
                                      REStats = RE_COVERAGE, GeneStats = GENE_COVERAGE,
                                      Seed = metadata(object)$Seed)
            
            return(object_trim)
          })

#' @rdname REMProduct-class
setMethod("rempAggregate", signature(object = "REMProduct"), 
          function(object, NCpG = 2) {
            REtype = object@REMPInfo[["REtype"]]
            method <- object@REMPInfo[["predictModel"]]
            
            if(grepl("aggregated", REtype))
            {
              message("The results are already aggregated. No changes made.")
              return(object) 
            }
            
            # If it is random forest prediction, then it has to be trimmed 
            # Otherwise it is fine.
            if(grepl("Random Forest", method))
            {
              if(!grepl("trimmed", method))
              {
                stop("Please first trim the results (using function 'rempTrim()').") 
              }
            }

            #The NCpG has to be >=2 
            if(NCpG < 2)
            {
              stop("The minimum number of CpGs must be >=2.")
            }
            
            M <- rempM(object)
            QC <- rempQC(object)
            RE_annotation <- rempAnnot(object)
            regionCode <-  metadata(object)$regionCode
            mcols(RE_annotation) <- cbind(mcols(RE_annotation), regionCode)
            GR <- rowRanges(object)
            
            ## Aggregate
            M_df <- DataFrame(M, Index = GR$RE.Index)
            M_df_MEAN <- aggregate(.~Index, M_df, function(x) mean(x, na.rm = TRUE))
            RE_annotation <- RE_annotation[match(M_df_MEAN$Index, RE_annotation$Index)]
            
            ## Remove RE with less than NCpG predicted
            count <- as.numeric(table(M_df$Index))
            M_df_MEAN <- M_df_MEAN[count>=NCpG,]
            
            M <- as.matrix(M_df_MEAN[,-1, drop = FALSE])

            ## Only works for random forest
            if(grepl("Random Forest", method))
            {
              QC_df <- DataFrame(QC, Index = GR$RE.Index)
              QC_df_MEAN <- aggregate(.~Index, QC_df, function(x) mean(x, na.rm = TRUE))
              QC_df_MEAN <- QC_df_MEAN[count>=NCpG,]
              QC <- as.matrix(QC_df_MEAN[,-1, drop = FALSE])
            } else {
              QC = NULL
            }
            
            RE_annotation <- RE_annotation[count>=NCpG]
            message(sum(count<NCpG), " ", REtype, " that have less than ", NCpG, " CpGs predicted are not aggretated.")
            
            RE_annotation_name <- colnames(mcols(RE_annotation))
            regionCode <- mcols(RE_annotation)[remp_options(".default.genomicRegionColNames")]
            RE_annotation <- RE_annotation[, RE_annotation_name[!RE_annotation_name %in% 
                                                                  remp_options(".default.genomicRegionColNames")]]
            cpgRanges <- granges(RE_annotation) # cpgRanges is now identical to RE ranges
            cpgRanges$RE.Index <- RE_annotation$Index

            refgene_main <- metadata(object)$refGene
            RE_CpG_ILMN <-  metadata(object)$RECpG
            
            ## Add rownames (RE index)
            RE_index <- as.character(cpgRanges$RE.Index)
            rownames(M) <- RE_index
            if(grepl("Random Forest", method))
            {
              rownames(QC) <- RE_index
            }

            ## Updated RE coverage
            RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN, 
                                             REtype, indent = "    ", FALSE)
            
            # Updated Gene coverage
            GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main, 
                                                 REtype, indent = "    ", FALSE)
            
            ## Update object
            object_aggregated <- REMProduct(REtype = paste0(REtype, " (aggregated by mean: min # of CpGs: ", NCpG, ")"), 
                                            platform = object@REMPInfo[["platform"]], 
                                            win = object@REMPInfo[["win"]],
                                            predictModel = method,  
                                            QCModel = object@REMPInfo[["QCModel"]], 
                                            rempM = M, rempQC = QC,
                                            cpgRanges = cpgRanges, sampleInfo = colData(object),
                                            REannotation = RE_annotation, 
                                            RECpG = RE_CpG_ILMN,
                                            regionCode = regionCode,
                                            refGene = metadata(object)$refGene,
                                            varImp = rempImp(object), 
                                            REStats = RE_COVERAGE, GeneStats = GENE_COVERAGE,
                                            Seed = metadata(object)$Seed)
            return(object_aggregated)
          })
