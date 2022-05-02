#' Main function for deconvolving signal using 3D contacts
#'
#' Optimizes the log-likelihood using BFGS and returns the deconvolved signal.
#'
#' @description This function fits a spatial beta-binomial model, estimates a
#' 2D smoothed probability function, and stores information for plotting results.
#' Uses a likelihood ratio test to test for a significant spatial fit over
#' optional baseline covariates such as cell type.
#'
#' @param y vector of observed counts
#'
#' @param D a symmetric matrix specifying the pairwise interactions. Should
#' be smoothed and normalized. Should be nxn where n is the length of y.
#' 
#' @param x (optional) vector of observed spatial positions. If not specified,
#' it is assumed that the counts y are evenly spaced.
#'
#' @param df integer, sets the number of degrees of freedom to use for the
#' smoothing spline. 
#' 
#' @param gamma (optional) weight for how much to penalize the off-diagonal
#' terms of D. default is 1, or no penalty. 
#' 
#' @return A list containing the
#' following output:
#' \itemize{ \item{\code{y}}{ a vector of the predicted deconvolved counts }
#' \item{\code{a}}{ a vector of the coefficient estimates}}.
#'
#'
#'
#' @export

decon <- function(y, D, x=NULL, df, gamma=1) {
  if (is.null(x)) {
    x <- seq(1,length(y))
  }
  B <- constructB(1:nrow(subchip), df=df)
  ILD <- (1-gamma)*diag(1,nrow=nrow(D))+gamma*D
  res <- optim(par=rep(0,df), loglik, y=y,D=ILD, B=B, method='BFGS')
  I <- solve(numDeriv::hessian(loglik, res$par, y=y, D=ILD, B=B))
  Sigma <- B%*%I%*%t(B)
}