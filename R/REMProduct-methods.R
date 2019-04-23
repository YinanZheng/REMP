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
  return(list(
    RE_Statistics = metadata(object)$REStats,
    Gene_Statistics = metadata(object)$GeneStats
  ))
})
# rempStats(object)


#### Utilities
#' @rdname REMProduct-class
setMethod(
  "plot", signature(x = "REMProduct", y = "missing"),
  function(x, type = c("individual", "overall"), ...) {
    type <- match.arg(type)
    default_main <- ""
    default_xlab <- "Methylation value (beta)"
    default_xlim <- c(0, 1)
    dotdotdot <- list(...)
    if (!hasArg("main")) {
      dotdotdot$main <- default_main
    }
    if (!hasArg("xlab")) {
      dotdotdot$xlab <- default_xlab
    }
    if (!hasArg("xlim")) {
      dotdotdot$xlim <- default_xlim
    }

    Bval <- rempB(x)

    if (type == "individual") {
      do.call(plot, c(list(x = density(na.omit(Bval[, 1]))), dotdotdot))
      if (ncol(Bval) >= 2) {
        for (i in seq(2, ncol(Bval)))
        {
          do.call(lines, c(list(x = density(na.omit(Bval[, i]))), dotdotdot))
        }
      }
    }

    if (type == "overall") {
      Bval.mean <- rowMeans(Bval, na.rm = TRUE)
      do.call(plot, c(list(x = density(na.omit(Bval.mean))), dotdotdot))
    }
  }
)
# plot(object)
# plot(object, type = "overall")

#' @rdname REMProduct-class
setMethod("details", signature(object = "REMProduct"), function(object) {
  REtype <- object@REMPInfo[["REtype"]]
  notAggregated <- TRUE
  if (grepl("aggregated", REtype)) notAggregated <- FALSE

  method <- object@REMPInfo[["predictModel"]]
  notProfiled <- TRUE
  notTrimmed <- TRUE
  notNaive <- TRUE
  if (grepl("Profiled", method)) notProfiled <- FALSE
  if (grepl("trimmed", method)) notTrimmed <- FALSE
  if (grepl("Naive", method)) notNaive <- FALSE

  .showREMPinfo(object)
  cat("\n")

  .showCpGcountbyChr(object)
  cat("\n")

  if (notProfiled & notAggregated & notTrimmed & notNaive) {
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
setMethod(
  "decodeAnnot", signature(object = "REMProduct"),
  function(object, type = c("symbol", "entrez"), ncore = 1, BPPARAM = NULL) {
    .isREMProductOrStop(object)

    be <- getBackend(ncore, BPPARAM, TRUE)

    skip <- TRUE
    type <- match.arg(type)
    message("Decoding ", object@REMPInfo[["REtype"]], " annotation to ", type, " ...")

    annot <- metadata(object)$REannotation
    code <- metadata(object)$regionCode
    refgene <- as.data.frame(metadata(object)$refGene)
    containSymbol <- any(grepl("symbol", colnames(mcols(annot))))
    containEntrez <- any(grepl("entrez", colnames(mcols(annot))))

    if (type == "symbol" & !containSymbol) {
      skip <- FALSE
      decodeFun <- .decodeSymbol
    }
    if (type == "entrez" & !containEntrez) {
      skip <- FALSE
      decodeFun <- .decodeEntrez
    }

    if (!skip) {
      BiocParallel::bpstart(be)
      .bploadLibraryQuiet("IRanges", be)
      # message("Pacakge loaded!")
      for (region in remp_options(".default.genomicRegionColNames")) {
        InRegion.original <- code[, region]
        ind.missing <- which(is.na(InRegion.original))
        InRegion <- na.omit(InRegion.original)
        res <- strsplit(InRegion, "[|]")
        res <- lapply(res, as.numeric)
        res <- IRanges::NumericList(res)
        res <- unlist(bplapply(res, decodeFun, refgene = refgene, BPPARAM = be))
        InRegion.original[-ind.missing] <- res
        annot$newregionEncode <- InRegion.original
        annot <- .changeColNames(
          annot, "newregionEncode",
          paste(region, type, sep = ".")
        )
      }
      metadata(object)$REannotation <- annot
      BiocParallel::bpstop(be)
    }
    return(object)
  }
)

#' @rdname REMProduct-class
setMethod(
  "rempTrim", signature(object = "REMProduct"),
  function(object, threshold = 1.7, missingRate = 0.2) {
    method <- object@REMPInfo[["predictModel"]]

    if (grepl("trimmed", method)) {
      previous_thres <- gsub("[\\(\\)]", "", regmatches(method, gregexpr("\\(.*?\\)", method))[[1]])
      if (threshold >= as.numeric(previous_thres)) {
        message("More stringent or equal threshold (", previous_thres, ") has been applied. No changes made.")
        return(object)
      } else {
        object@REMPInfo[["predictModel"]] <- method <- "Random Forest"
      }
    }

    if (!grepl("Random Forest", method)) {
      message("rempTrim() is only applicable to prediction using Random Forest model. No changes made.")
      return(object)
    }

    REtype <- object@REMPInfo[["REtype"]]
    # beta <- as.matrix(rempB(object))
    M <- as.matrix(rempM(object))
    QC <- as.matrix(rempQC(object))
    badInd <- QC > threshold
    QC[badInd] <- NA
    # beta[badInd] <- NA
    M[badInd] <- NA

    removeCpG <- rowMeans(badInd) > missingRate

    RE_annotation <- rempAnnot(object)
    RE_CpG_ILMN <- metadata(object)$RECpG
    regionCode <- metadata(object)$regionCode
    refgene_main <- metadata(object)$refGene

    # beta <- beta[!removeCpG,,drop = FALSE]
    M <- M[!removeCpG, , drop = FALSE]
    QC <- QC[!removeCpG, , drop = FALSE]
    cpgRanges <- rowRanges(object)[!removeCpG, ]

    mcols(RE_annotation) <- cbind(mcols(RE_annotation), regionCode)
    RE_annotation <- subsetByOverlaps(RE_annotation, cpgRanges) # RE_annotation is updated here
    RE_annotation_name <- colnames(mcols(RE_annotation))
    regionCode <- mcols(RE_annotation)[remp_options(".default.genomicRegionColNames")]
    RE_annotation <- RE_annotation[, RE_annotation_name[!RE_annotation_name %in%
      remp_options(".default.genomicRegionColNames")]]

    ## Updated RE coverage
    RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN,
      REtype,
      indent = "    ", FALSE
    )

    # Updated Gene coverage
    GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main,
      REtype,
      indent = "    ", FALSE
    )


    ## Update object
    object_trim <- REMProduct(
      REtype = object@REMPInfo[["REtype"]],
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
      Seed = metadata(object)$Seed
    )

    return(object_trim)
  }
)

#' @rdname REMProduct-class
setMethod(
  "rempAggregate", signature(object = "REMProduct"),
  function(object, NCpG = 2, ncore = 1, BPPARAM = NULL) {
    be <- getBackend(ncore, BPPARAM, TRUE)

    REtype <- object@REMPInfo[["REtype"]]
    method <- object@REMPInfo[["predictModel"]]

    if (grepl("aggregated", REtype)) {
      message("The results are already aggregated. No changes made.")
      return(object)
    }

    # If it is random forest prediction, then it has to be trimmed
    # Otherwise it is fine.
    if (grepl("Random Forest", method)) {
      if (!grepl("trimmed", method)) {
        stop("Please first trim the results (using function 'rempTrim()').")
      }
    }

    # The NCpG has to be >=2
    if (NCpG < 2) {
      stop("The minimum number of CpGs must be >=2.")
    }

    M <- rempM(object)
    QC <- rempQC(object)
    RE_annotation <- rempAnnot(object)
    regionCode <- metadata(object)$regionCode
    mcols(RE_annotation) <- cbind(mcols(RE_annotation), regionCode)
    GR <- rowRanges(object)

    ## Remove RE with less than NCpG predicted
    M_df <- data.frame(M,
      Index = as.character(GR$RE.Index),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    RE_annotation <- RE_annotation[match(
      unique(M_df$Index),
      runValue(RE_annotation$Index)
    )]
    # identical(runValue(RE_annotation$Index), unique(as.character(M_df$Index)))
    select_RE <- runValue(RE_annotation$Index)
    count_table <- table(M_df$Index)
    count_table <- count_table[select_RE]
    # identical(names(count_table), as.character(RE_annotation$Index))
    count <- as.numeric(count_table)
    RE_annotation <- RE_annotation[count >= NCpG]
    select_RE <- select_RE[count >= NCpG]
    # identical(select_RE, as.character(RE_annotation$Index))
    M_df <- M_df[M_df$Index %in% select_RE, ]

    ## Aggregate predicted M values
    if (ncore > 1) {
      BiocParallel::bpstart(be)
      ITER <- .iblkrow_dup(M_df, chunks = ncore, index_col_name = "Index")
      M_df_MEAN <- BiocParallel::bpiterate(ITER, .aggregateREMP,
        index_col_name = "Index",
        BPPARAM = be,
        REDUCE = rbind,
        reduce.in.order = TRUE
      )
      BiocParallel::bpstop(be)
    } else {
      M_df_MEAN <- aggregate(. ~ Index, M_df, mean, na.rm = TRUE, na.action = na.pass)
      M_df_MEAN <- M_df_MEAN[match(select_RE, M_df_MEAN$Index), ]
    }
    # identical(M_df_MEAN_par, M_df_MEAN)
    # identical(runValue(RE_annotation$Index), as.character(M_df_MEAN$Index))
    M <- as.matrix(M_df_MEAN[, -1, drop = FALSE])

    ## Aggregate QC only works for random forest
    if (grepl("Random Forest", method)) {
      QC_df <- data.frame(QC,
        Index = GR$RE.Index, check.names = FALSE,
        stringsAsFactors = FALSE
      )
      QC_df <- QC_df[QC_df$Index %in% select_RE, ]

      if (ncore > 1) {
        BiocParallel::bpstart(be)
        ITER <- .iblkrow_dup(QC_df, chunks = ncore, index_col_name = "Index")
        QC_df_MEAN <- BiocParallel::bpiterate(ITER, .aggregateREMP,
          index_col_name = "Index",
          BPPARAM = be,
          REDUCE = rbind,
          reduce.in.order = TRUE
        )
        BiocParallel::bpstop(be)
      } else {
        QC_df_MEAN <- aggregate(. ~ Index, QC_df, mean, na.rm = TRUE, na.action = na.pass)
        QC_df_MEAN <- QC_df_MEAN[match(select_RE, QC_df_MEAN$Index), ]
      }
      # identical(QC_df_MEAN_par, QC_df_MEAN)
      QC <- as.matrix(QC_df_MEAN[, -1, drop = FALSE])
    } else {
      QC <- NULL
    }
    # identical(QC_df_MEAN$Index, M_df_MEAN$Index)
    # identical(as.character(RE_annotation$Index), as.character(M_df_MEAN$Index))

    message(sum(count < NCpG), " ", REtype, " that have less than ", NCpG, " CpGs predicted are not aggretated.")

    ## Wrapping up
    RE_annotation_name <- colnames(mcols(RE_annotation))
    regionCode <- mcols(RE_annotation)[remp_options(".default.genomicRegionColNames")]
    RE_annotation <- RE_annotation[, RE_annotation_name[!RE_annotation_name %in%
      remp_options(".default.genomicRegionColNames")]]
    cpgRanges <- granges(RE_annotation) # cpgRanges is now identical to RE ranges
    cpgRanges$RE.Index <- RE_annotation$Index

    refgene_main <- metadata(object)$refGene
    RE_CpG_ILMN <- metadata(object)$RECpG

    ## Add rownames (RE index)
    RE_index <- runValue(cpgRanges$RE.Index)
    rownames(M) <- RE_index
    if (grepl("Random Forest", method)) {
      rownames(QC) <- RE_index
    }

    ## Updated RE coverage
    RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN,
      REtype,
      indent = "    ", FALSE
    )

    # Updated Gene coverage
    GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main,
      REtype,
      indent = "    ", FALSE
    )

    ## Update object
    object_aggregated <- REMProduct(
      REtype = paste0(REtype, " (aggregated by mean: min # of CpGs: ", NCpG, ")"),
      platform = object@REMPInfo[["platform"]],
      win = object@REMPInfo[["win"]],
      predictModel = method,
      QCModel = object@REMPInfo[["QCModel"]],
      rempM = M, rempQC = QC,
      cpgRanges = cpgRanges, sampleInfo = colData(object),
      REannotation = RE_annotation,
      RECpG = RE_CpG_ILMN,
      regionCode = regionCode,
      refGene = refgene_main,
      varImp = rempImp(object),
      REStats = RE_COVERAGE, GeneStats = GENE_COVERAGE,
      Seed = metadata(object)$Seed
    )
    return(object_aggregated)
  }
)

#' @rdname REMProduct-class
setMethod(
  "rempCombine", signature(object1 = "REMProduct", object2 = "REMProduct"),
  function(object1, object2) {
    REtype1 <- object1@REMPInfo[["REtype"]]
    method1 <- object1@REMPInfo[["predictModel"]]
    platform1 <- object1@REMPInfo[["platform"]]
    win1 <- object1@REMPInfo[["win"]]
    seed1 <- metadata(object1)$Seed
    cnames1 <- colnames(object1)

    REtype2 <- object2@REMPInfo[["REtype"]]
    method2 <- object2@REMPInfo[["predictModel"]]
    platform2 <- object2@REMPInfo[["platform"]]
    win2 <- object2@REMPInfo[["win"]]
    seed2 <- metadata(object2)$Seed
    cnames2 <- colnames(object2)

    duplicated_samples <- intersect(cnames1, cnames2)

    if (!identical(REtype1, REtype2)) stop(paste0("You cannot combine ", REtype1, " methylation data with ", REtype2, " methylation data."))
    if (!identical(method1, method2)) stop(paste0("You cannot combine ", method1, " prediction with ", method2, " prediction."))
    if (!identical(platform1, platform2)) stop(paste0("You cannot combine ", platform1, " array with ", platform2, " array."))
    if (!identical(win1, win2)) stop(paste0("You cannot combine prediction with flanking window size = ", win1, " with size = ", win2, "."))
    if (!identical(seed1, seed2)) stop(paste0("You cannot combine prediction with using seed = ", win1, " with seed = ", win2, "."))
    if (length(duplicated_samples) > 0) stop(paste0("Duplicated samples detected: ", paste0(duplicated_samples, collapse = ", ")))

    M1 <- as.matrix(rempM(object1))
    M2 <- as.matrix(rempM(object2))
    nrow1 <- nrow(M1)
    nrow2 <- nrow(M2)
    if (nrow1 != nrow2) stop("Number of rows does mot match!")

    B1 <- as.matrix(rempB(object1))
    B2 <- as.matrix(rempB(object2))

    QC1 <- as.matrix(rempQC(object1))
    QC2 <- as.matrix(rempQC(object2))

    QCModel <- object1@REMPInfo[["QCModel"]]

    cpgranges1 <- rowRanges(object1)
    sampleinfo1 <- colData(object1)
    reannot1 <- rempAnnot(object1)
    regionCode1 <- metadata(object1)$regionCode
    mcols(reannot1) <- cbind(mcols(reannot1), regionCode1)
    RE_CpG_ILMN1 <- metadata(object1)$RECpG
    varimp1 <- rempImp(object1)

    cpgranges2 <- rowRanges(object2)
    sampleinfo2 <- colData(object2)
    reannot2 <- rempAnnot(object2)
    regionCode2 <- metadata(object2)$regionCode
    mcols(reannot2) <- cbind(mcols(reannot2), regionCode2)
    RE_CpG_ILMN2 <- metadata(object2)$RECpG
    varimp2 <- rempImp(object2)

    ## Combine:
    M <- cbind(M1, M2)
    B <- cbind(B1, B2)
    QC <- cbind(QC1, QC2)
    if (nrow(QC) == 0) QC <- NULL

    cpgRanges <- cpgranges1
    sampleinfo <- rbind(sampleinfo1, sampleinfo2)
    varimp <- DataFrame(varimp1, varimp2, check.names = FALSE)

    RE_annotation <- unique(c(reannot1, reannot2))
    RE_CpG_ILMN <- unique(c(RE_CpG_ILMN1, RE_CpG_ILMN2))

    ## Wrapping up
    RE_annotation_name <- colnames(mcols(RE_annotation))
    regionCode <- mcols(RE_annotation)[remp_options(".default.genomicRegionColNames")]
    RE_annotation <- RE_annotation[, RE_annotation_name[!RE_annotation_name %in%
      remp_options(".default.genomicRegionColNames")]]

    refgene_main <- metadata(object1)$refGene

    ## Updated RE coverage
    RE_COVERAGE <- .coverageStats_RE(RE_annotation, regionCode, cpgRanges, RE_CpG_ILMN,
      REtype1,
      indent = "    ", FALSE
    )

    # Updated Gene coverage
    GENE_COVERAGE <- .coverageStats_GENE(regionCode, refgene_main,
      REtype1,
      indent = "    ", FALSE
    )


    ## Update object
    object_combined <- REMProduct(
      REtype = REtype1,
      platform = platform1,
      win = win1,
      predictModel = method1,
      QCModel = object1@REMPInfo[["QCModel"]],
      rempM = M, rempB = B, rempQC = QC,
      cpgRanges = cpgRanges, sampleInfo = sampleinfo,
      REannotation = RE_annotation,
      RECpG = RE_CpG_ILMN,
      regionCode = regionCode,
      refGene = refgene_main,
      varImp = varimp,
      REStats = RE_COVERAGE, GeneStats = GENE_COVERAGE,
      Seed = seed1
    )
    return(object_combined)
  }
)
