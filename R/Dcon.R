#' constructor of \code{\linkS4class{Dcon}} object
#' @param counts An integer vector of counts for each bin
#' @param bins An integer vector of bin positions, usually chromosome coordinates
#' in base pairs.
#' @param D a sparse matrix representing the contacts between bins (e.g. Hi-C
#' matrix)
#'
#' Counts should be untransformed count-level data
#'
#' @return Returns a \code{\linkS4class{Dcon}} object containing the counts,
#' bins, and 
#' @export
Dcon <- function(counts, cell_types, nUMI = NULL, require_int = TRUE, n_max_cells = 10000) {
  counts <- check_counts(counts, 'Reference', require_2d = T, require_int = require_int)
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'Reference', require_2d = T, require_int = require_int)
  }
  check_cell_types(cell_types)
  barcodes <- intersect(intersect(names(nUMI), names(cell_types)), colnames(counts))
  if(length(barcodes) == 0)
    stop('Reference: cell_types, counts, and nUMI do not share any barcode names. Please ensure that names(cell_types) matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),length(cell_types),dim(counts)[2]))
    warning('Reference: some barcodes in nUMI, cell_types, or counts were not mutually shared. Such barcodes were removed.')
  if(sum(nUMI[barcodes] != colSums(counts[,barcodes])) > 0)
    warning('Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this is intended, there is no problem.')
  missing_cell_types <- names(which(table(cell_types[barcodes]) == 0))
  if(length(missing_cell_types) > 0)
    warning(paste('Reference: missing cell types with no occurences: ',paste(missing_cell_types,collapse=', ')))
  reference <- new("Reference", cell_types = cell_types[barcodes], counts = counts[,barcodes], nUMI = nUMI[barcodes])
  cur_count <- max(table(reference@cell_types))
  if(cur_count > n_max_cells) {
    warning(paste0('Reference: number of cells per cell type is ', cur_count, ', larger than maximum allowable of ', n_max_cells,
                   '. Downsampling number of cells to: ', n_max_cells))
    reference <- create_downsampled_data(reference, n_samples = n_max_cells)
  }
  return(reference)
}