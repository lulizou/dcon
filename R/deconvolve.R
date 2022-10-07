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
          def = function(object, df=NULL, gamma=0.5, pad=10, save_mat=F) {
            if (is.null(df)) {
              df <- pad*2+1
            }
            coords <- object@coords
            loops <- object@loops
            hic <- object@hic
            if (length(hic)==1) {
              hic <- hic[[1]]
            }
            hichip_mats <- object@mat
            if (is.null(names(hichip_mats))) {
              stop('hichip matrices must be input as a named list')
            }
            results <- before <- after <- NULL
            if (save_mat) {
              new_mat <- list(before = list(), after = list())
            }
            if (nrow(loops)==0) {
              message('No loops found - assuming that each matrix in mat is 
                      a loop window')
              if (length(hichip_mats)>1) {
                pb <- txtProgressBar(min = 1, max = length(hichip_mats),
                                     initial = 1, style = 3)
              }
              for (i in 1:length(hichip_mats)) {
                if (i>1) {
                  setTxtProgressBar(pb,i)
                }
                m <- hichip_mats[[i]]
                r1 <- rowSums(m)
                r2 <- colSums(m)
                hic1 <- hic[[1]]
                hic2 <- hic[[2]]
                B <- construct_basis(1:length(r1), df=df)
                fit_r1 <- try(fit_decon(r1, hic1), silent=T)
                fit_r2 <- try(fit_decon(r2, hic2), silent=T)
                if (is(fit_r1, 'try-error') | is(fit_r2, 'try-error')) {
                  next
                }
                a1 <- fit_r1$a
                a2 <- fit_r2$a
                ff1 <- (B%*%a1)-log(sum(exp(B%*%a1)))
                ff2 <- (B%*%a2)-log(sum(exp(B%*%a2)))
                res <- exp(matrix(data = rep(ff1, df), nrow=df, ncol=df) + 
                             t(matrix(data = rep(ff2, df), nrow = df, ncol = df)))*sum(m)
                res[res<1e-2] <- 0
                if (save_mat) {
                  new_mat$before[[names(hichip_mats)[i]]] <- m
                  new_mat$after[[names(hichip_mats)[i]]] <- res
                }
              }
            } else {
              message('Found loops - assuming that each matrix in mat is 
                      a hichip matrix')
              before <- after <- matrix(nrow = nrow(loops), ncol = length(hichip_mats))
              colnames(before) <- colnames(after) <- names(hichip_mats)
              pb <- txtProgressBar(min = 1, max = nrow(loops), initial = 1, style=3) 
              chrom_start <- object@coords$id[1]
              for (i in 1:nrow(loops)) {
                setTxtProgressBar(pb,i)
                a1 <- loops$bin1[i]
                a2 <- loops$bin2[i]
                a1_idx <- seq(a1-pad,a1+pad)
                a2_idx <- seq(a2-pad,a2+pad)
                hic1_idx <- a1_idx-chrom_start+1
                hic2_idx <- a2_idx-chrom_start+1
                if (any(c(hic1_idx, hic2_idx)>min(dim(hic))) | any(c(hic1_idx,hic2_idx)<0)) {
                  next
                }
                hic1 <- normalize_hic(hic[hic1_idx,hic1_idx], gamma=gamma)
                hic2 <- normalize_hic(hic[hic2_idx,hic2_idx], gamma=gamma)
                myres <- list()
                for (j in 1:length(hichip_mats)) {
                  if (any(c(a1_idx,a2_idx)>min(dim(hichip_mats[[j]])))| any(c(a1_idx,a2_idx)<0)) {
                    next
                  }
                  m <- hichip_mats[[j]][a1_idx,a2_idx]
                  before[i,j] <- sum(m[pad:(pad+2),pad:(pad+2)])
                  if (nnzero(m)> 0) {
                    r1 <- rowSums(m)
                    r2 <- colSums(m)
                    B <- construct_basis(1:length(r1), df=df)
                    fit_r1 <- try(fit_decon(r1, hic1), silent=T)
                    fit_r2 <- try(fit_decon(r2, hic2), silent=T)
                    if (is(fit_r1, 'try-error') | is(fit_r2, 'try-error')) {
                      next
                    }
                    a1 <- fit_r1$a
                    a2 <- fit_r2$a
                    ff1 <- (B%*%a1)-log(sum(exp(B%*%a1)))
                    ff2 <- (B%*%a2)-log(sum(exp(B%*%a2)))
                    res <- exp(matrix(data = rep(ff1, df), nrow=df, ncol=df) + 
                                 t(matrix(data = rep(ff2, df), nrow = df, ncol = df)))*sum(m)
                    res[res<1e-2] <- 0
                    myres[[names(hichip_mats)[j]]] <- res
                    after[i,j] <- sum(res[pad:(pad+2),pad:(pad+2)])
                    if (save_mat) {
                      new_mat$before[[paste0('loop_',i,'_',names(hichip_mats)[j])]] <- m
                      new_mat$after[[paste0('loop_',i,'_',names(hichip_mats)[j])]] <- res
                    }
                  }
                }
                if (length(myres)>2) {
                  diffs <- as.data.frame(rbind(
                    get_pairwise_diff(hichip_mats, a1_idx, a2_idx, label='before'),
                    get_pairwise_diff(myres, 1:(pad*2+1), 1:(pad*2+1), label='after')
                  ))
                  diffs$loop <- i
                  if (is.null(results)) {
                    results <- diffs
                  } else {
                    results <- rbind(results, diffs)
                  }
                }
              }
              close(pb)
              object@loops <- loops
            }
            if (save_mat) {
              object@mat <- new_mat
            } else {
              object@mat <- list()
            }
            object@results <- list(df = results, before = before, after = after)
            object@is_deconvolved <- TRUE
            return(object)
          })