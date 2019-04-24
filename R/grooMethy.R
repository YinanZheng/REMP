#' @title Groom methylation data to fix potential data issues
#'
#' @description
#' \code{grooMethy} is used to automatically detect and fix data issues including zero beta
#' value, missing value, and infinite value.
#'
#' @param methyDat A \code{\link{RatioSet}}, \code{\link{GenomicRatioSet}}, \code{\link{DataFrame}},
#' data.table, data.frame, or matrix of Illumina BeadChip methylation data (450k or EPIC array).
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
#' For methylation data in beta value, if zero value exists, the logit transformation
#' from beta to M value will produce negative infinite value. Therefore, zero beta value
#' will be replaced with the smallest non-zero beta value found in the dataset. \code{grooMethy}
#' can also handle missing value (i.e. \code{NA} or \code{NaN}) using KNN imputation (see
#' \code{\link{impute.knn}}). The infinite value will be also treated as missing value for imputation.
#' If the original dataset is in beta value, \code{grooMethy} will first transform it to M value
#' before imputation is carried out. If the imputed value is out of the original range (which is possible when
#' \code{imputebyrow = FALSE}), mean value will be used instead. Warning: imputed
#' values for multimodal distributed CpGs (across samples) may not be correct. Please check package \code{ENmix} to
#' identify the CpGs with multimodal distribution. Please note that \code{grooMethy} is
#' also embedded in \code{\link{remp}} so the user can run \code{\link{remp}} directly without
#' explicitly running \code{grooMethy}.
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
grooMethy <- function(methyDat, impute = TRUE, imputebyrow = TRUE, mapGenome = FALSE, verbose = FALSE) {
  currenT <- Sys.time()
  methyDat_work <- methyDat

  if (is(methyDat_work, "RatioSet") || is(methyDat_work, "GenomicRatioSet")) {
    methyDat_work <- minfi::getBeta(methyDat_work)
  }

  methyDat_work <- .methyMatrix(methyDat_work)

  type <- .guessBetaorM(methyDat_work)
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
  if (arrayType == "UNKNOWN") {
    stop("Unknown methylation array type.")
  }

  if (verbose) {
    message(
      "Illumina ", arrayType, " Methylation data in ", type,
      " value detected."
    )
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

    ## If any beta value = 0 is found
    if (2 %in% dc$code) {
      if (verbose) {
        message(
          "    A total of ", sum(methyDat_work == 0, na.rm = TRUE),
          " zero beta values are found."
        )
      }
      smallBeta <- min(methyDat_work[methyDat_work > 0], na.rm = TRUE)
      if (verbose) {
        message(
          "    Fixing zero beta values with the smallest beta value = ",
          smallBeta, " ..."
        )
      }
      methyDat_work[methyDat_work == 0] <- smallBeta
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
      message("    Converting beta value to M value ...")
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
        message("    Fixing ", nrow(idx), " imputed probes that are out of the original data range ....")
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

  rset <- minfi::RatioSet(Beta = betadata, M = mdata, annotation = annotationInfo)
  if (mapGenome) {
    rset <- minfi::mapToGenome(rset)
  }

  if (verbose) message("Methylation data grooming is completed.", .timeTrace(currenT)$t_text)

  return(rset)
} ## End of grooMethy

## Internal functions
.methyMatrix <- function(methyDat) {
  methyDat <- data.frame(methyDat, check.names = FALSE)
  probeNameIndicator <- which(vapply(methyDat, class, character(1)) %in% c("factor", "character"))
  if (length(probeNameIndicator) > 1) {
    stop("There are more than one character columns! Please only keep one column or just use row names to indicate Illumina probe names (i.e. cg00000029).")
  }

  containILMN <- "cg" %in% unique(substring(methyDat[seq_len(10), probeNameIndicator], 1, 2))
  containRownames <- "cg" %in% unique(substring(rownames(methyDat)[seq_len(10)], 1, 2))
  if (!containILMN & !containRownames) {
    stop("Cannot find a column or row names that indicates Illumina probe names (i.e. cg00000029). Please fix it.")
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
  return(methyDat.matrix)
}

.dataCheck <- function(methyDat, type) {
  errorString <- NULL
  code <- NULL
  if (type %in% c("[Genomic]RatioSet", "beta") & any(methyDat < 0 | methyDat >
    1, na.rm = TRUE)) {
    errorString <- c(errorString, "out-of-range")
    code <- c(code, 1)
  }
  if (type %in% c("[Genomic]RatioSet", "beta") & any(methyDat == 0, na.rm = TRUE)) {
    errorString <- c(errorString, "0")
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
