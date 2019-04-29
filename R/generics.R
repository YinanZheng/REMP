### REMParcel-class

#' @rdname REMParcel-class
#' @exportMethod getParcelInfo
setGeneric("getParcelInfo", function(object) standardGeneric("getParcelInfo"))

#' @rdname REMParcel-class
#' @exportMethod getRefGene
setGeneric("getRefGene", function(object) standardGeneric("getRefGene"))

#' @rdname REMParcel-class
#' @exportMethod getRE
setGeneric("getRE", function(object) standardGeneric("getRE"))

#' @rdname REMParcel-class
#' @exportMethod getRECpG
setGeneric("getRECpG", function(object) standardGeneric("getRECpG"))

#' @rdname REMParcel-class
#' @exportMethod getILMN
setGeneric("getILMN", function(object, ...) standardGeneric("getILMN"))

#' @rdname REMParcel-class
#' @exportMethod saveParcel
setGeneric("saveParcel", function(object, ...) standardGeneric("saveParcel"))


### REMProduct-class

#' @rdname REMProduct-class
#' @exportMethod rempM
setGeneric("rempM", function(object) standardGeneric("rempM"))

#' @rdname REMProduct-class
#' @exportMethod rempB
setGeneric("rempB", function(object) standardGeneric("rempB"))

#' @rdname REMProduct-class
#' @exportMethod rempQC
setGeneric("rempQC", function(object) standardGeneric("rempQC"))

#' @rdname REMProduct-class
#' @exportMethod rempAnnot
setGeneric("rempAnnot", function(object) standardGeneric("rempAnnot"))

#' @rdname REMProduct-class
#' @exportMethod rempImp
setGeneric("rempImp", function(object) standardGeneric("rempImp"))

#' @rdname REMProduct-class
#' @exportMethod rempStats
setGeneric("rempStats", function(object) standardGeneric("rempStats"))

#' @rdname REMProduct-class
#' @exportMethod remplot
setGeneric("remplot", function(object, ...) standardGeneric("remplot"))

#' @rdname REMProduct-class
#' @exportMethod details
setGeneric("details", function(object) standardGeneric("details"))

#' @rdname REMProduct-class
#' @exportMethod decodeAnnot
setGeneric("decodeAnnot", function(object, ...) standardGeneric("decodeAnnot"))

#' @rdname REMProduct-class
#' @exportMethod rempTrim
setGeneric("rempTrim", function(object, ...) standardGeneric("rempTrim"))

#' @rdname REMProduct-class
#' @exportMethod rempAggregate
setGeneric("rempAggregate", function(object, ...) standardGeneric("rempAggregate"))

#' @rdname REMProduct-class
#' @exportMethod rempCombine
setGeneric("rempCombine", function(object1, object2) standardGeneric("rempCombine"))
