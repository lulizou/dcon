#' S4 constructors for Dcon object 
#' 
#' 

Dcon <- setClass(
  'Dcon',
  
  slots = c(
    mat = 'list',
    coords = 'GRanges',
    anchors = 'GRanges',
    anchors_hic = 'list'
  ),
  
  prototype=list(
    mat = list(),
    coords = GRanges(),
    anchors = GRanges(),
    anchors_hic = list()
  ),
  
  validity = function(object) {
    if (is.null(object@coords$id)) {
      return('coords must have 1 metadata column for the id')
    }
    if (length(anchors) != length(anchors_hic)) {
      return('anchors and anchors_hic must be the same length')
    }
    return(TRUE)
  }
)

setGeneric(
  name = 'deconvolve',
  def = function(object) {
    standardGeneric('deconvolve')
  }
)

setMethod(
  f = 'deconvolve',
  def = function(object) {
    
  }
)




