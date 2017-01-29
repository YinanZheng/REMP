### REMParcel-class

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
#' @exportMethod annotation
setGeneric("annotation", function(object) standardGeneric("annotation"))

#' @rdname REMProduct-class
#' @exportMethod imp
setGeneric("imp", function(object) standardGeneric("imp"))

#' @rdname REMProduct-class
#' @exportMethod stats
setGeneric("stats", function(object) standardGeneric("stats"))

#' @rdname REMProduct-class
#' @exportMethod details
setGeneric("details", function(object) standardGeneric("details"))

#' @rdname REMProduct-class
#' @exportMethod decodeAnnot
setGeneric("decodeAnnot", function(object, ...) standardGeneric("decodeAnnot"))

#' @rdname REMProduct-class
#' @exportMethod trim
setGeneric("trim", function(object, ...) standardGeneric("trim"))
