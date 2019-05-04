#' @title Groom methylation data to fix potential data issues
#'
#' @description
#' \code{grooMethy} is used to automatically detect and fix data issues including zero beta
#' value, missing value, and infinite value.
#'
#' @param methyDat A \code{\link{RatioSet}}, \code{\link{GenomicRatioSet}}, \code{\link{DataFrame}},
#' \code{data.table}, \code{data.frame}, or \code{matrix} of Illumina BeadChip methylation data
#' (450k or EPIC array) or Illumina methylation sequencing data.
#' @param Seq.GR A \code{\link{GRanges}} object containing genomic locations of the CpGs profiled by sequencing
#' platforms. This parameter should not be \code{NULL} if the input methylation data \code{methyDat} are
#' obtained by sequencing. Note that the genomic location must be in hg19 build. See details.
#' @param impute If \code{TRUE}, K-Nearest Neighbouring imputation will be applied to fill
#' the missing values. Default = \code{TRUE}. See Details.
#' @param imputebyrow If \code{TRUE}, missing values will be imputed using similar values in row
#' (i.e., across samples); if \code{FALSE}, missing values will be imputed using similar values
#' in column (i.e., across CpGs). Default is \code{TRUE}.
#' @param mapGenome Logical parameter. If \code{TRUE}, function will return a \code{\link{GenomicRatioSet}}
#' object instead of a \code{\link{RatioSet}}.
#' @param verbose Logical parameter. Should the function be verbose?
#'
#' @details
#' For methylation data in beta value, if zero/one value exists, the logit transformation
#' from beta to M value will produce infinite value. Therefore, zero/one beta value
#' will be replaced with the smallest non-zero beta/largest non-one beta value found in the dataset.
#' \code{grooMethy} can also handle missing value (i.e. \code{NA} or \code{NaN}) using KNN imputation (see
#' \code{\link{impute.knn}}). The infinite value will be also treated as missing value for imputation.
#' If the original dataset is in beta value, \code{grooMethy} will first transform it to M value
#' before imputation is carried out. If the imputed value is out of the original range (which is possible when
#' \code{imputebyrow = FALSE}), mean value will be used instead. Warning: imputed
#' values for multimodal distributed CpGs (across samples) may not be correct. Please check package \code{ENmix} to
#' identify the CpGs with multimodal distribution. Please note that \code{grooMethy} is
#' also embedded in \code{\link{remp}} so the user can run \code{\link{remp}} directly without
#' explicitly running \code{grooMethy}. For sequencing methylation data, please specify the genomic location of CpGs
#' in a \code{GenomicRanges} object and specify it in \code{Seq.GR}. For an example of \code{Seq.GR}, Please
#' run \code{minfi::getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)} (the row names of the CpGs in \code{Seq.GR}
#' can be \code{NULL}).
#'
#' @return A \code{\link{RatioSet}} or \code{\link{GenomicRatioSet}} containing beta value and
#' M value of the methylation data.
#'
#' @examples
#' # Get GM12878 methylation data (450k array)
#' if (!exists("GM12878_450k")) GM12878_450k <- getGM12878("450k")
#' GM12878_450k <- grooMethy(GM12878_450k, verbose = TRUE)
#' 
#' # Also works if data input is a matrix
#' grooMethy(minfi::getBeta(GM12878_450k), verbose = TRUE)
#' @export
grooMethy <- function(methyDat, Seq.GR = NULL, impute = TRUE, imputebyrow = TRUE, mapGenome = FALSE, verbose = FALSE) {
  currenT <- Sys.time()

  if (is.null(methyDat)) stop("Methylation dataset (methyDat) is missing.")
  if (!is.null(Seq.GR) & !is(Seq.GR, "GRanges")) stop("Seq.GR must be a GenomicRanges object.")

  methyDat_work <- methyDat

  if (is(methyDat_work, "RatioSet") || is(methyDat_work, "GenomicRatioSet")) {
    if(minfi::preprocessMethod(methyDat_work) == "grooMethy(REMP)") {
      if (verbose) message("Methylation data are already groomed.")
      return(methyDat_work)
    }  
    methyDat_work <- minfi::getBeta(methyDat_work)
  }

  methyDat_work <- .methyMatrix(methyDat_work, Seq.GR)

  type <- .guessDataType(methyDat_work)
  arrayType <- .guessArrayType(methyDat_work)

  if (arrayType == "27k") {
    stop("Illumina 27k array is not supported.")
  }
  if (arrayType == "450k") {
    annotationInfo <- c(array = "IlluminaHumanMethylation450k", annotation = remp_options(".default.450k.annotation"))
  }
  if (arrayType == "EPIC") {
    annotationInfo <- c(array = "IlluminaHumanMethylationEPIC", annotation = remp_options(".default.epic.annotation"))
  }
  if (arrayType == "Sequencing" | !is.null(Seq.GR)) {
    annotationInfo <- c(array = "IlluminaHumanMethylationSequencing", annotation = "Custom")
  }

  if (verbose) {
    message(
      "Illumina ", arrayType, " Methylation data in ", type,
      " value detected."
    )
  }

  if (type == "percentage") {
    methyDat_work <- methyDat_work / 100
    if (verbose) message("Percentage data have been divided by 100 and scaled to range 0-1.")
  }

  dc <- .dataCheck(methyDat_work, type)

  if (is.null(dc$code)) {
    if (verbose) {
      message("    No issue found in the methylation dataset.")
    }
    impute <- FALSE
  } else {
    ## If beta and any beta value is out of range
    if (1 %in% dc$code) {
      stop("Negative beta or beta > 1 values detected! Please check your data and data preprocessing.")
    }

    ## If any beta value = 0 or 1 is found
    if (2 %in% dc$code) {
      nzero <- sum(methyDat_work == 0, na.rm = TRUE)
      none <- sum(methyDat_work == 1, na.rm = TRUE)

      if (nzero > 0) {
        if (verbose) message("    A total of ", nzero, " zero beta values are found.")
        smallBeta <- min(methyDat_work[methyDat_work > 0], na.rm = TRUE)
        if (verbose) message("    Fixing 'zero' beta values with the smallest non-0 beta value = ", smallBeta, " ...")
        methyDat_work[methyDat_work == 0] <- smallBeta
      }

      if (none > 0) {
        if (verbose) message("    A total of ", none, " one beta values are found.")
        bigBeta <- max(methyDat_work[methyDat_work < 1], na.rm = TRUE)
        if (verbose) message("    Fixing 'one' beta values with the largest non-1 beta value = ", bigBeta, " ...")
        methyDat_work[methyDat_work == 1] <- bigBeta
      }
    }

    ## If any beta/M value is NaN or infinite
    if (any(c(3, 4, 5) %in% dc$code)) {
      methyDat_work[is.nan(methyDat_work)] <- NA
      methyDat_work[is.infinite(methyDat_work)] <- NA

      if (verbose) {
        message("    A total of ", sum(is.na(methyDat_work)), " NA/NaN/Inf values are found.")
      }
    } else {
      impute <- FALSE
    }
  }

  ## Prepare both beta and M value using corrected data
  if (type == "M") {
    if (verbose) {
      message("    Converting M value to beta value ...")
    }
    betadata <- .toBeta(methyDat_work)
    mdata <- methyDat_work
  } else {
    if (verbose) {
      message("    Converting beta/percentage value to M value ...")
    }
    betadata <- methyDat_work
    mdata <- .toM(methyDat_work)
  }

  ## Impute the data!
  if (impute) {
    if (verbose) {
      message("    Imputing missing values using KNN method ...")
    }
    if (ncol(mdata) == 1) {
      stop("KNN-imputation cannot be applied to single sample.")
    }

    if (imputebyrow) resu <- impute::impute.knn(t(mdata)) else resu <- impute::impute.knn(mdata)

    ## Error checking and imperfect fix (use mean to replace the
    ## out-of-range imputation)
    rg <- apply(mdata, 1, function(x) range(x, na.rm = TRUE)) # Calculate the range of the data before imputation

    if (imputebyrow) imputed_mdata <- t(resu$data) else imputed_mdata <- resu$data

    ## Check if any imputed data is out of the original range
    idx <- which(imputed_mdata < rg[1, ] | imputed_mdata > rg[2, ],
      arr.ind = TRUE
    )

    ## If so then use mean to fix them
    if (nrow(idx) > 0) {
      if (verbose) {
        message("    Fixing ", nrow(idx), " imputed probes that are out of the original data range...")
      }
      m <- apply(mdata[idx[, 1], ], 1, function(x) mean(x, na.rm = TRUE)) ## calculate the mean of the original data
      for (i in seq_len(nrow(idx))) {
        imputed_mdata[idx[i, 1], idx[i, 2]] <- m[i]
      }
    } else {
      if (verbose) {
        message("    All imputed probes are within the original data range.")
      }
    }

    ## Update with imputed data
    mdata <- imputed_mdata
    betadata <- .toBeta(mdata)
  }

  rset <- minfi::RatioSet(Beta = betadata, M = mdata, 
                          annotation = annotationInfo, 
                          preprocessMethod = "grooMethy(REMP)")
  if (arrayType != "Sequencing" & mapGenome) {
    rset <- minfi::mapToGenome(rset)
  }

  if (verbose) message("Methylation data grooming is completed.", .timeTrace(currenT)$t_text)

  return(rset)
} ## End of grooMethy

## Internal functions
.methyMatrix <- function(methyDat, Seq.GR = NULL) {
  methyDat <- data.frame(methyDat, check.names = FALSE)

  probeNameIndicator <- which(vapply(methyDat, class, character(1)) %in% c("factor", "character"))

  if (!is.null(Seq.GR)) {
    if (length(probeNameIndicator) > 0) {
      stop("Factor or Character columns detected! The input methylation data from sequencing platform should all be numeric.")
    }
    methyDat.matrix <- as.matrix(methyDat)
    rownames(methyDat.matrix) <- paste0(seqnames(Seq.GR), ":", start(Seq.GR))
  } else {
    if (length(probeNameIndicator) > 1) {
      stop(paste(
        "For array methylation data, please only keep one column or just use row names to indicate Illumina probe names (i.e. cg00000029).",
        "For sequencing methylation data, the parameter Seq.GR cannot be missing. Please provide it."
      ))
    }

    containILMN <- "cg" %in% unique(substring(methyDat[seq_len(10), probeNameIndicator], 1, 2))
    containRownames <- "cg" %in% unique(substring(rownames(methyDat)[seq_len(10)], 1, 2))
    if (!containILMN & !containRownames) {
      stop(paste(
        "For array methylation data, a column or row names that indicates Illumina probe names (i.e. cg00000029) is missing.",
        "Please fix it. For sequencing methylation data, the parameter Seq.GR cannot be missing. Please provide it."
      ))
    }

    if (containILMN) {
      methyDat.matrix <- matrix(methyDat[, -probeNameIndicator],
        ncol = ncol(methyDat) - 1
      )
      rownames(methyDat.matrix) <- methyDat[, probeNameIndicator]
      colnames(methyDat.matrix) <- colnames(methyDat)[-probeNameIndicator]
    } else if (containRownames) {
      methyDat.matrix <- as.matrix(methyDat)
    }
  }
  return(methyDat.matrix)
}

.dataCheck <- function(methyDat, type) {
  errorString <- NULL
  code <- NULL
  if (type %in% c("[Genomic]RatioSet", "beta", "percentage") & any(methyDat < 0 | methyDat >
    1, na.rm = TRUE)) {
    errorString <- c(errorString, "out-of-range")
    code <- c(code, 1)
  }
  if (type %in% c("[Genomic]RatioSet", "beta", "percentage") & any(methyDat == 0 | methyDat == 1, na.rm = TRUE)) {
    errorString <- c(errorString, "0/1")
    code <- c(code, 2)
  }
  if (any(is.nan(methyDat))) {
    errorString <- c(errorString, "NaN")
    code <- c(code, 3)
  }
  if (any(is.na(methyDat))) {
    errorString <- c(errorString, "NA")
    code <- c(code, 4)
  }
  if (any(is.infinite(methyDat))) {
    errorString <- c(errorString, "Inf")
    code <- c(code, 5)
  }
  return(list(code = code, errorString = errorString))
}
