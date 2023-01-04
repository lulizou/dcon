#' Normalize Hi-C or DNA-DNA contact matrix before deconvolution
#' 
#' Column-normalizes and optionally smooths
#' 
#' @param M the matrix of raw counts
#' @param gamma the amount of penalty on the off-diagonal between 0-1. Default is
#' 1 (no penalty) 
#' @param threshold the minimum number of counts in a column.
#' If a column has less than this number, it does not participate in 
#' deconvolution. Default is 5 (but usually should go higher)
#' @param smooth boolean indicating whether or not to smooth the contact 
#' matrix. Smoothing is accomplished with `spatstat::blur`. The operations
#' are moderately time intensive and this probably needs to be sped up at 
#' some point, so default is `FALSE`.
#' @param sigma if smooth is `TRUE`, the sd of the Gaussian kernel smooth. Default
#' is 1. May need to tune this if experimenting with different sized M.
#' 
#' @return a list containing `a`, the estimated coefficients, `B`, the
#' basis functions, `est` the estimated log-scale signal and `var` the variance
#' associated with the estimates.
#' 
#' @importFrom spatstat.core blur
#' @importFrom spatstat.geom im
#' 
#' @export

normalize_hic <- function(M, gamma=1, threshold = 5, smooth = F, sigma = 1) {
  # threshold exists so that low count columns don't have large influence
  M <- as.matrix(M)
  if (any(is.na(M))) {
    stop('there are NA values in your matrix!')
  }
  if (!any(diag(M)>0)) {
    stop('there are no diagonal entries on the matrix that are non-zero')
  }
  thresh_idx <- which(colSums(M)<=threshold)
  if (length(thresh_idx)>0) {
    for (id in thresh_idx) {
      M[,id] <- 0 # 0 out the column
      M[id,id] <- mean(diag(M)) # set to the mean of the diagonal
    }
  }
  if (smooth) {
    M <- as.matrix(blur(im(M), sigma = sigma))
  }
  Dsmoothnorm <- sweep(M,2,colSums(M),'/')
  ILD <- (1-gamma)*diag(1,nrow=nrow(Dsmoothnorm))+gamma*Dsmoothnorm
  return(ILD)
}