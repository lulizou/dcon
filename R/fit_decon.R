#' Fit the 1D deconvolution model
#' 
#' Fits the 1D version of the deconvolution model, where a 1D
#' vector of counts gets convolved by a 2D matrix
#' 
#' @param reads a vector of counts 
#' @param D a (normalized) matrix representing the DNA-DNA contacts. See
#' `normalize_hic` to normalize a matrix of counts. Must be symmetric with
#' number of rows and columns equal to the length of reads.
#' @param df the degrees of freedom for the Gaussian basis functions. The default
#' is the length of reads.
#' @param offset a vector of offset counts the same length as reads.
#' 
#' @return a list containing `a`, the estimated coefficients, `B`, the
#' basis functions, `est` the estimated log-scale signal and `var` the variance
#' associated with the estimates.
#' 
#' @export

fit_decon <- function(reads, D, df=NULL, offset=NULL) {
  if (is.null(df)) {
    df <- length(reads)
  }
  B <- construct_basis(1:length(reads), df=df)
  if (is.null(offset)) {
    starting_a <- rep(0, df)
  } else {
    starting_a <- rep(-mean(log(offset)), df)
  }
  if (any(offset==0)) {
    stop('some values of offset are 0; consider omitting or adding 1')
  }
  if (!is.null(offset)) {
    D <- offset*D
  }
  o <- optim(par=starting_a, loglik, y=reads, D=D, B=B, method='BFGS')
  a <- o$par
  # get hessian
  H <- hess_a(y=reads, D=D, B=B, a=a)
  Sigma <- solve(psd(H))
  BSigmaB <- B%*%Sigma%*%t(B)
  Bavars <- diag(BSigmaB)
  return(list(a=a, B=B, est=as.numeric(B%*%a), var = Bavars))
}