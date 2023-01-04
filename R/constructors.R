#' S4 constructors for Dcon object 
#' 
#' @import GenomicRanges
#' @import Matrix
#' 
#' @export

Dcon <- setClass(
  
  'Dcon',
  
  slots = c(
    mat = 'list',
    coords = 'GRanges',
    loops = 'data.frame',
    hic = 'Matrix',
    is_deconvolved = 'logical',
    results = 'list'
  ),
  
  prototype=list(
    mat = list(),
    coords = GRanges(),
    loops = data.frame(),
    hic = Matrix(),
    is_deconvolved = F,
    results = list()
  ),
  
  validity = function(object) {
    if (is.null(names(object@mat))) {
      return('hichip matrices must be input as a named list')
    }
    if (is.null(coords$id)) {
      
    }
    return(TRUE)
  }
)


setGeneric(
  name = 'deconvolve',
  def = function(object, ...) {
    standardGeneric('deconvolve')
  }
)

setMethod(
  f = 'show',
  signature = 'Dcon',
  definition = function(object) {
    cat('Dcon object\n')
    cat(paste0('Number of matrices: ', length(object@mat)))
    if (length(object@coords)>0) {
      cat(paste0('\nResolution: ', object@coords[1]@ranges@width-1, 'bp'))
    }
    cat(paste0('\nNumber of loops: ', nrow(object@loops)))
    cat(paste0('\nHiC dimensions: ', paste0(dim(object@hic),collapse='x')))
    if (object@is_deconvolved) {
      cat('\nDeconvolved')
    } else {
      cat('\nNot deconvolved')
    }
  }
)






