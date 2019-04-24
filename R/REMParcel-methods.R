
# REMParcel-methods

setMethod("show", signature(object = "REMParcel"), function(object) {
  .showREMParceInfo(object)
})
# remparcel

#' @rdname REMParcel-class
setMethod("saveParcel", signature(object = "REMParcel"), function(object, work.dir = tempdir(), verbose = FALSE, ...) {
  arrayType <- object@REMParcelInfo[["platform"]]
  REType <- object@REMParcelInfo[["REtype"]]
  win <- object@REMParcelInfo[["max.win"]]

  subDirName <- paste0("REMP.data.", arrayType)
  work.dir <- .forwardSlashPath(work.dir)
  data.dir <- .forwardSlashPath(file.path(work.dir, subDirName))

  if (dir.create(data.dir, showWarnings = FALSE)) {
    if (verbose) {
      message(
        "A new directory ", paste0("'", subDirName, "'"),
        " has been created under the directory:\n", data.dir
      )
    }
  } else {
    if (verbose) {
      message(
        "Directory ", paste0("'", subDirName, "'"),
        " already exists under the directory:\n", data.dir
      )
    }
  }

  path_to_parcel <- file.path(
    data.dir,
    paste0(
      "REMParcel_",
      REType, "_",
      arrayType, "_",
      win, ".rds"
    )
  )

  saveRDS(object, file = path_to_parcel, ...)

  message("REMParcel has been saved under the directory:\n", data.dir)
})

#' @rdname REMParcel-class
setMethod("getParcelInfo", signature(object = "REMParcel"), function(object) {
  return(object@REMParcelInfo)
})

#' @rdname REMParcel-class
setMethod("getRefGene", signature(object = "REMParcel"), function(object) {
  return(object@RefGene)
})

#' @rdname REMParcel-class
setMethod("getRE", signature(object = "REMParcel"), function(object) {
  return(object@RE)
})

#' @rdname REMParcel-class
setMethod("getRECpG", signature(object = "REMParcel"), function(object) {
  return(object@RECpG)
})

#' @rdname REMParcel-class
setMethod("getILMN", signature(object = "REMParcel"), function(object, REonly = FALSE) {
  ilmn <- object@ILMN
  if (REonly) {
    return(ilmn[!is.na(ilmn$RE.Index)])
  } else {
    return(ilmn)
  }
})
