#'
#' Run deconvolution on a set of loops specified in a `Dcon` object
#'
#' @description Description
#'
#' @param object the `Dcon` object containing the HiChIP matrix or matrices,
#' coordinates, loops, and HiC/DNA-DNA contacts
#' @param df the degrees of freedom for the basis functions; default will
#' set it to the maximum i.e. the length of the window
#' @param gamma between 0 and 1, controls whether or not the off-diagonal is
#' penalized (0) or not (1). default is 0.5
#' @param pad the amount of padding to add around the loop center to create
#' the window for deconvolution
#' @param save_mat boolean indicating whether or not to save each deconvolved
#' loop window in the results field. Default is F.
#' @param store_sun boolean indicating whether or not to store the sum of all
#' the before/after matrices. Used to compute the average deconvolved loop
#' window after deconvolution.
#' 
#' @return A `Dcon` object which now also has a results field containing the
#' results of the deconvolution. 
#'
#'
#'
#' @export
setMethod(f = 'deconvolve',
          signature = 'Dcon',
          def = function(object, df=NULL, gamma=0.5, pad=10, save_mat=F,
                         store_sum = F) {
            win_len <- pad*2+1
            if (is.null(df)) {
              df <- win_len
            }
            coords <- object@coords
            loops <- object@loops
            hic <- object@hic
            hichip_mats <- object@mat
            before <- after <- matrix(nrow = nrow(loops), ncol = length(hichip_mats))
            colnames(before) <- colnames(after) <- names(hichip_mats)
            results <- NULL
            if (save_mat) {
              new_mat <- list(before = list(), after = list())
            }
            if (store_sum) {
              sum_before <- matrix(data = NA, nrow = win_len, ncol = win_len)
              sum_after <- matrix(data = NA, nrow = win_len, ncol = win_len)
            }
            ndconvolved <- 0
            if (nrow(loops)>1) {
              pb <- txtProgressBar(min = 1, max = nrow(loops), initial = 1, style=3) 
            }
            for (i in 1:nrow(loops)) {
              if (nrow(loops)>1) {
                setTxtProgressBar(pb,i)
              }
              a1 <- loops$bin1[i]
              a2 <- loops$bin2[i]
              a1_idx <- seq(a1-pad,a1+pad)
              a2_idx <- seq(a2-pad,a2+pad)
              if (any(a1_idx<0) | any(a2_idx<0)) {
                next
              }
              hic1 <- hic[a1_idx,a1_idx]
              hic2 <- hic[a2_idx,a2_idx]
              if (!any(diag(hic1)>0) | !any(diag(hic2)>0)) {
                next
              }
              hic1 <- normalize_hic(hic1, gamma=gamma)
              hic2 <- normalize_hic(hic2, gamma=gamma)
              myres <- list()
              for (j in 1:length(hichip_mats)) {
                m <- hichip_mats[[j]][a1_idx,a2_idx]
                before[i,j] <- sum(m[pad:(pad+2),pad:(pad+2)])
                if (nnzero(m)> 0) {
                  r1 <- rowSums(m)
                  r2 <- colSums(m)
                  B <- construct_basis(1:length(r1), df=df)
                  fit_r1 <- fit_decon(r1, hic1, df=df)
                  fit_r2 <- fit_decon(r2, hic2, df=df)
                  a1 <- fit_r1$a
                  a2 <- fit_r2$a
                  ff1 <- (B%*%a1)-log(sum(exp(B%*%a1)))
                  ff2 <- (B%*%a2)-log(sum(exp(B%*%a2)))
                  res <- exp(matrix(data = rep(ff1, length(ff2)), nrow=length(ff1), ncol=length(ff2)) + 
                               t(matrix(data = rep(ff2, length(ff1)), nrow = length(ff2), ncol = length(ff1))))*sum(m)
                  res[res<1e-2] <- 0
                  myres[[names(hichip_mats)[j]]] <- res
                  after[i,j] <- sum(res[pad:(pad+2),pad:(pad+2)])
                  if (save_mat) {
                    new_mat$before[[paste0('loop_',i,'_',names(hichip_mats)[j])]] <- m
                    new_mat$after[[paste0('loop_',i,'_',names(hichip_mats)[j])]] <- res
                  }
                  if (store_sum) {
                    if (any(is.na(sum_before))) {
                      sum_before <- m
                      sum_after <- res
                    } else {
                      sum_before <- sum_before+m
                      sum_after <- sum_after+res
                    }
                    ndconvolved <- ndconvolved + 1
                  }
                }
              }
              if (length(myres)>=2) {
                # get the pairwise differences for the whole matrix
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
            if (nrow(loops)>1) {
              close(pb)
            }
            if (save_mat) {
              object@mat <- new_mat
            } else if (store_sum) {
              object@mat <- list(sum_before = sum_before, sum_after = sum_after, 
                                 n = ndconvolved)
            } else {
              object@mat <- list()
            }
            object@results <- list(df = results, before = before, after = after)
            object@is_deconvolved <- TRUE
            return(object)
          })