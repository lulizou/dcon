#'
#' Optimizes the log-likelihood using BFGS and returns the deconvolved signal.
#'
#' @description Description
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
setMethod(f = 'deconvolve',
          signature = 'Dcon',
          def = function(object, df=21, gamma=0.5, pad=10, add_center=1) {
            center_idx <- median(seq(1,2*pad+1))
            center_idx <- seq(center_idx-add_center,center_idx+add_center)
            loops <- object@loops
            hichip_mats <- object@mat
            after <- matrix(data=NA, nrow = nrow(loops), ncol = length(hichip_mats))
            colnames(after) <- names(hichip_mats)
            before <- after
            pb <- txtProgressBar(min = 1, max = nrow(loops), initial = 1, style=3) 
            chrom_start <- object@coords$id[1]
            for (i in 1:nrow(loops)) {
              setTxtProgressBar(pb,i)
              a1 <- loops$bin1[i]
              a2 <- loops$bin2[i]
              a1_idx <- seq(a1-pad,a1+pad)
              a2_idx <- seq(a2-pad,a2+pad)
              hic1_idx <- a1_idx-cs+1
              hic2_idx <- a2_idx-cs+1
              hic1 <- normalize_hic(hic[hic1_idx,hic1_idx], gamma=gamma)
              hic2 <- normalize_hic(hic[hic2_idx,hic2_idx], gamma=gamma)
              for (j in 1:length(hichip_mats)) {
                m <- hichip_mats[[j]][a1_idx,a2_idx]
                before[i,j] <- sum(m[center_idx,center_idx])
                if (nnzero(m)> 0) {
                  r1 <- rowSums(m)
                  r2 <- colSums(m)
                  B <- construct_basis(1:length(r1), df=df)
                  fit_r1 <- try(fit_decon(r1, hic1))
                  fit_r2 <- try(fit_decon(r2, hic2))
                  if (is(fit_r1, 'try-error') | is(fit_r2, 'try-error')) {
                    next
                  }
                  a1 <- fit_r1$a
                  a2 <- fit_r2$a
                  ff1 <- (B%*%a1)-log(sum(exp(B%*%a1)))
                  ff2 <- (B%*%a2)-log(sum(exp(B%*%a2)))
                  res <- exp(matrix(data = rep(ff1, df), nrow=df, ncol=df) + 
                               t(matrix(data = rep(ff2, df), nrow = df, ncol = df)))*sum(m)
                  after[i,j] <- sum(res[center_idx,center_idx])
                  res[res<1e-2] <- 0
                  # hichip_mats[[j]][a1_idx,a2_idx] <- res
                }
              }
            }
            close(pb)
            object@mat <- hichip_mats
            object@results <- list(before=before, after=after)
            object@is_deconvolved <- TRUE
            return(object)
          })