# accessory functions for importing data


#' Get pairwise L1 and L2 norm differences from a list of matrices
#' and the indices at those matrices to calculate the differences on
get_pairwise_diff <- function(matrix_list, idx1, idx2, label='') {
  diffs <- data.frame(t(combn(names(matrix_list),2)))
  colnames(diffs) <- c('m1','m2')
  diffs$label <- label
  diffs$diff <- NA
  for (i in 1:nrow(diffs)) {
    diffs$diff[i] <- sum(abs(matrix_list[[diffs$m1[i]]][idx1, idx2]-matrix_list[[diffs$m2[i]]][idx1,idx2]))
  }
  return(diffs)
}

#' Load in anchors
#' If there are multiple anchor files, can decide whether to use all unique
#' entries or only take ones in both (intersection). Default is to take the
#' intersection.
#' @param ... one or more file names to input. 
#' @param type the type of the anchor file. Supports fithichip
#' @param use whether or not to use the intersection or union if there are 
#' multiple input files
load_anchors <- function(..., coords, type='fithichip', use='union',
                         resolution = NULL, chrom = NULL) {
  if (!use%in%c('union','intersection')) {
    stop('use parameter must be either union or intersection')
  }
  if (!type%in%c('fithichip','hichipper')) {
    stop('use parameter must be either fithichip or hichipper')
  }
  if (type == 'hichipper' & is.null(resolution)) {
    stop('if loading in hichipper, must specify resolution of loops e.g. 2500')
  }
  input_files <- list(...)
  dfs <- list()
  bin_pairs <- NULL
  for (i in 1:length(input_files)) {
    d <- fread(input_files[[i]])
    if (type == 'hichipper') {
      colnames(d) <- c('chr1', 's1', 'e1', 'chr2', 's2', 'e2', 'idk1', 'count')
      d$source <- input_files[[i]]
      d$s1 <- round(d$s1/resolution)*resolution
      d$s2 <- round(d$s2/resolution)*resolution
      d$e1 <- d$s1+resolution
      d$e2 <- d$s2+resolution
    }
    if (!is.null(chrom)) {
      d <- d %>% filter(chr1==chrom)
    }
    d$bin1 <- left_join(
      d[,c('chr1','s1','e1')] %>% 
        dplyr::rename(seqnames=chr1, start=s1, end=e1),
      as.data.frame(coords) %>% 
        dplyr::select(seqnames, start, end, id), 
      by = c("seqnames", "start", "end")) %>% 
      pull(id)
    d$bin2 <- left_join(d[,c('chr2','s2','e2')] %>% 
                          dplyr::rename(seqnames=chr2, start=s2, end=e2),
                        as.data.frame(coords) %>% 
                          dplyr::select(seqnames, start, end, id),
                        by = c("seqnames", "start", "end")) %>% 
      pull(id)
    d$bin_pair <- paste0(d$bin1,',',d$bin2)
    dfs[[i]] <- d
    if (is.null(bin_pairs)) {
      bin_pairs <- d$bin_pair
    } else {
      if (use == 'intersection') {
        bin_pairs <- bin_pairs[bin_pairs %in% dfs[[i]]$bin_pair]
      } 
    }
  }
  if (use == 'intersection') {
    dfs[[1]] <- dfs[[1]] %>% dplyr::filter(bin_pair %in% bin_pairs)
    return(dfs[[1]])
  } else if (use == 'union') {
    big_df <- dfs[[1]]
    for (i in 2:length(dfs)) {
      big_df <- bind_rows(big_df, dfs[[i]])
    }
    return(big_df %>% distinct())
  }
  
}


#' Load in the hic-pro coordinates file
#' aka the _abs.bed file with the correct resolution bins
load_coords <- function(filename, type='hic-pro') {
  if (type == 'hic-pro') {
    coords <- fread(filename, col.names = c('chrom','start','end','id'))
  } else if (type == 'windows') {
    coords <- fread(filename, col.names=c('chrom','start','end'))
    coords$id <- seq(1,nrow(coords))
  }
  coords <- makeGRangesFromDataFrame(coords, keep.extra.columns = T)
  return(coords)
}


#' Load in the hichip matrix as a a sparseMatrix
#' type specifies what type the file is. Currently only supports the hic-pro
#' .matrix format. Should be made from allValidPairs. 
#' @param filename specifies the name of the file to load in
#' @param coords specifies the total dimensions of the matrix, created from bins
#' using load_coords.
#' @param chrom specifies which chromosome to subset to based on coords. Default
#' is NULL so it will return the whole thing.
#' @param type specifies the file type
#' @param diagonal specifies which entries to 0 out from the diagonal.
load_hichip <- function(filename, coords, chrom=NULL,
                        type='hic-pro', diagonal=3) {
  d <- fread(filename)
  d <- d %>% filter(abs(V2-V1)>diagonal)
  d <- sparseMatrix(i = d$V1, j = d$V2, x = d$V3, dims=c(length(coords),length(coords)))
  if (!is.null(chrom)) {
    coords_sub <- coords[seqnames(coords)==chrom,]
    d <- d[coords_sub$id, coords_sub$id]
  }
  return(d)
}





# functions for specifically hichip

# symmetrically pad on either side of each grange
pad_windows <- function(granges, window_size=5e4, bin_size=2500) {
  x <- (window_size/2)-floor(width(ranges(granges))/2)
  starts <- start(ranges(granges))-x
  ends <- end(ranges(granges))+x
  short_idx <- which((ends-starts)!=window_size)
  ends[short_idx] <- ends[short_idx]+1 
  # add one to the starts that are divisible so it doesn't pick up the previous
  starts_0 <- (starts%%bin_size)==0
  starts[starts_0] <- starts[starts_0]+1
  start(ranges(granges)) <- starts
  end(ranges(granges)) <- ends
  granges <- trim(granges)
  return(granges)
}
normalize_hic <- function(M, gamma=1, threshold = 5) {
  # threshold exists so that low count columns won't get deconvolved 
  M <- as.matrix(M)
  thresh_idx <- which(colSums(M)<=threshold)
  if (length(thresh_idx)>0) {
    for (id in thresh_idx) {
      M[,id] <- 0 # 0 out the column
      M[id,id] <- 1 # make the diagonal = 1
    }
  }
  Dsmoothnorm <- sweep(M,2,colSums(M),'/')
  ILD <- (1-gamma)*diag(1,nrow=nrow(Dsmoothnorm))+gamma*Dsmoothnorm
  return(ILD)
}

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

decon_to_2d <- function(decon1, decon2, mat) {
  new_mat <- (decon1$est/sum(decon1$est))%*%(t(decon2$est)/sum(decon2$est))
  new_mat <- new_mat*sum(mat)
  return(new_mat)
}



# Get gaussian basis functions over an interval (vector)
construct_basis <- function(interval, df = NULL) {
  if (is.null(df)) {
    df <- length(interval)
  }
  if (df < 3) {
    stop('Interval and df must be >= 5')
  }
  range <- diff(range(interval))
  kn <- seq(min(interval) - (range/(df-3))*3, max(interval) + (range/(df))*3, 
            by = range/(df-3))
  myu <- kn[3:(df+2)]
  h <- diff(kn,lag = 2)/3
  B <- matrix(0,length(interval),(df))
  for (j in 1:df){
    B[,j] <- exp(-0.5*(interval-myu[j])^2/(h[1]^2))
  }
  return(B)
}

#' Simulate the log(rate) using a linear combination of Gaussian basis functions
#' assuming that the coefficients are normally distributed with specified
#' mean and variance
simulate_log_rate <- function(len, df, alpha_mean = -1, alpha_sd = 2,
                              npeaks = NULL) {
  B <- construct_basis(1:len, df = df)
  if (is.null(npeaks)) {
    a <- rnorm(n = df, mean = alpha_mean, sd = alpha_sd)
  } else {
    a <- rep(-2, len)
    peak_locations <- 1+c(1:npeaks)*floor(len/(npeaks+1))
    a[peak_locations] <- 5
  }
  return(list(B = B, a = a))
}

#' Simulate a convolving matrix that resembles a typical DNA-DNA contact matrix
#' decay is the ratio between each consecutive diagonal
#' add_noise is the rate parameter for adding Poisson noise
#' if tads > 1, adds in equal-length tads by dividing len/tads
simulate_D <- function(len, decay = 1/3, gamma = 1, add_noise = NULL,
                       tads = 1, tad_strength = 0.2) {
  D <- diag(len)
  for (i in 1:len) {
    D[row(D) == (col(D)-i)] <- (1/i)*decay
    D[row(D) == (col(D)+i)] <- (1/i)*decay
  }
  if (tads > 1) {
    tad_size <- floor(len/tads)
    for (i in 1:tads) {
      p1 <- (i-1)*tad_size+1
      p2 <- i*tad_size
      D[p1:p2,p1:p2] <- D[p1:p2,p1:p2]+tad_strength
    }
  }
  if (!is.null(add_noise)) {
    if (!is.numeric(add_noise)) {
      stop('add_noise parameter must be either NULL or a numeric value
           corresponding to the rate parameter of the poisson noise')
    }
    D <- D*add_noise*5
    D <- D + matrix(rpois(len*len, add_noise), nrow = len, ncol = len)
  }
  
  return(normalize_hic(D, gamma = gamma))
}

#' Simulate the Poisson-distributed outcome y 
simulate_y <- function(len, df, D = NULL, alpha_mean = -1, alpha_sd = 2,
                       decay = 1/3, gamma = 1, npeaks = NULL) {
  log_rate_params <- simulate_log_rate(len, df, alpha_mean, alpha_sd, npeaks = npeaks)
  if (is.null(D)) {
    D <- simulate_D(len, decay, gamma)
  }
  B <- log_rate_params$B
  a <- log_rate_params$a
  y <- rpois(n = len, lambda = D%*%exp(B%*%a))
  return(list(y = y, D = D, B = B, a = a))
}


myquadprog <- function(a0, y, D, B, conv = 0.001) {
  lam_threshold = 1e-8; init_val <- -5; MIN_ITERATIONS <- 6
  beta_succ <- 1.1; beta_fail <- 0.5; gamma <- 0.1
  epsilon_2 <- 1e-5; delta <- 0.1 #trust region radius
  current_a <- a0
  prev_loglik <- loglik(current_a,y,D,B)
  diff <- 100
  while (diff >= conv) {
    prev_a <- current_a
    H_0 <- hess_a(y,D,B,current_a)
    grad_0 <- -deriv_a(y,D,B,current_a)
    H <- psd(H_0)
    #norm_factor <- norm(H,'2')
    #H <- H/norm_factor
    #grad <- grad_0/norm_factor
    #epsilon <- 1e-7
    #H <- H + epsilon*diag(length(grad))
    grad <- grad_0
    A <- diag(length(grad))
    solution <- quadprog::solve.QP(H,grad,A,current_a,meq=0)$solution
    predicted_decrease <- -(0.5*t(solution)%*%H_0%*%solution - sum(grad_0*solution))
    print(solution)
    current_loglik <- loglik(solution,y,D,B)
    true_decrease <- prev_loglik-current_loglik
    print(prev_loglik)
    print(current_loglik)
    print(true_decrease)
    print(gamma*predicted_decrease)
    if (true_decrease >= gamma*predicted_decrease) {
      # then we expand the trust region
      delta <- beta_succ*delta
      current_a <- solution
      prev_loglik <- current_loglik
    } else {
      # otherwise decrease the trust region
      delta <- min(1, beta_fail*delta)
    }
    diff <- norm(as.matrix(prev_a-current_a))
    print(prev_a)
    print(current_a)
    message(diff)
  }
  return(current_a)
}
psd <- function(H) {
  eig <- eigen(H); epsilon = 1e-3
  if(length(H) == 1)
    P <- eig$vectors %*% pmax(eig$values,epsilon) %*% t(eig$vectors)
  else
    P <- eig$vectors %*% diag(pmax(eig$values,epsilon)) %*% t(eig$vectors)
  return(P)
}
loglik <- function(a, y, D, B) {
  # get the negative log-likelihood
  eBa <- exp(B%*%a)
  DeBa <- D%*%eBa
  first_term <- sum(y*log(DeBa))
  second_term <- sum(DeBa)
  third_term <- sum(lfactorial(y))
  return(-first_term+second_term+third_term)
  
}
ll_negbin <- function(y, mu, theta) {
  return(lgamma(y+theta)-lgamma(theta)-lgamma(y+1)+theta*log(theta)
         -theta*log(mu+theta)+y*log(mu)-y*log(mu+theta))
}
loglik_negbin <- function(a,theta, y, D, B) {
  mu <- D%*%exp(B%*%a)
  return(-sum(ll_negbin(y, mu, theta)))
}
# loglik_negbin_optim <- function(atheta, y, D, B) {
#   n <- length(y)
#   npar <- ncol(B)
#   a <- atheta[1:npar]
#   theta <- atheta[length(atheta)]
#   eBa <- exp(B%*%a)
#   DeBa <- D%*%eBa
#   first_term <- sum(lgamma(y+theta))
#   second_term <- n*lgamma(theta)
#   third_term <- sum(lgamma(y+1))
#   fourth_term <- sum(theta*(log(theta)))
#   fifth_term <- sum(theta*log(DeBa))
#   sixth_term <- theta^2
#   seventh_term <- sum(y*log(DeBa+theta))
#   eighth_term <- sum(y*log(DeBa))
#   return(-first_term + second_term + third_term - fourth_term + fifth_term
#          - sixth_term + seventh_term - eighth_term)
# }
deriv_a_negbin <- function(y, D, B, a, theta) {
  eBa <- exp(B%*%a)
  DeBa <- D%*%eBa
  eBab <- sweep(B, 1, FUN='*', eBa)
  DeBab <- D%*%eBab
  first_term <- colSums(sweep(DeBab, 1, FUN='*', theta*(1/(DeBa+theta))))
  second_term <- colSums(sweep(DeBab, 1, FUN='*', y*(1/(DeBa))))
  third_term <- colSums(sweep(DeBab, 1, FUN='*', y*(1/(DeBa+theta))))
  return(-first_term+second_term-third_term)
}
d_theta_negbin <- function(y, mu, theta) {
  return(digamma(y+theta)-digamma(theta)+log(theta)+1-log(mu+theta)-theta/(mu+theta)
         -y/(mu+theta))
}
deriv_theta_negbin <- function(y, D, B, a, theta) {
  mu <- D%*%exp(B%*%a)
  return(sum(d_theta_negbin(y,mu,theta)))
}
grad_descent_negbin <- function(y, D, B, theta, a_init, step=1e-6, conv = 0.001,
                                maxiter=5e3, verbose=F) {
  current_a <- a_init
  prev_loglik <- loglik_negbin(current_a,theta,y,D,B)
  diff <- 100
  ni<-0
  while (diff >= conv) {
    ni<-ni+1
    prev_a <- current_a
    grad_0 <- deriv_a_negbin(y,D,B,current_a,theta)
    current_a <- prev_a + grad_0*step
    current_loglik <- loglik_negbin(current_a,theta,y,D,B)
    diff <- norm(as.matrix(prev_a-current_a))
    if (verbose) {
      message(diff)
    }
    if(ni==maxiter) {
      message(paste0('Reached maxiter ', maxiter))
      break
    }
  }
  return(current_a)
}
grad_descent_negbin_both <- function(y, D, B, theta, a_init, step=1e-3, conv = 0.001,
                                     maxiter=5e3, verbose=F) {
  current_a <- a_init
  current_theta <- theta
  prev_loglik <- loglik_negbin(current_a,current_theta,y,D,B)
  prev_diff <- 1e6
  diff <- 100
  ni<-0
  while (diff >= conv) {
    ni<-ni+1
    prev_a <- current_a
    prev_theta <- current_theta
    grad_0 <- deriv_a_negbin(y,D,B,current_a,current_theta)
    grad_theta <- deriv_theta_negbin(y,D,B,current_a,current_theta)
    current_a <- prev_a + grad_0*step
    current_theta <- prev_theta + grad_theta*step
    current_loglik <- loglik_negbin(current_a,current_theta,y,D,B)
    diff <- norm(as.matrix(prev_loglik-current_loglik))
    current_diff <- diff
    prev_loglik <- current_loglik
    if (verbose) {
      message(diff)
    }
    if(ni==maxiter) {
      message(paste0('Reached maxiter ', maxiter))
      break
    }
  }
  return(list(a=current_a,theta=current_theta))
}

partial_ak <- function(y, D, B, a, k) {
  # compute the partial derivative of the loglik w.r.t. a_k
  n <- length(y)
  deriv_result <- 0
  bk <- B[,k,drop=F]
  for (j in 1:n) {
    Dj <- D[,j,drop=F]
    term1 <- as.numeric(t(Dj)%*%exp(B%*%a))
    term2 <- as.numeric(t(Dj)%*%(exp(B%*%a)*bk))
    deriv_result <- deriv_result + y[j]*(1/term1)*term2 - term2
  }
  return(deriv_result)
}
deriv_a <- function(y, D, B, a) {
  eBa <- exp(B%*%a)
  DeBa <- D%*%eBa
  eBab <- sweep(B, 1, FUN='*', eBa)
  DeBab <- D%*%eBab
  first_term <- colSums(sweep(DeBab, 1, FUN='*', y*(1/DeBa)))
  second_term <- colSums(DeBab)
  return(second_term-first_term)
}
hess_a <- function(y, D, B, a) {
  eBa <- exp(B%*%a)
  DeBa <- D%*%eBa
  DeBa <- DeBa+1e-9
  eBab <- sweep(B, 1, FUN='*', eBa)
  DeBab <- D%*%eBab
  # fix j
  result <- matrix(NA, nrow=length(a), ncol=length(a))
  for (j in 1:length(a)) {
    eBab_j <- eBa*B[,j]
    # first term
    first_term <- t(sweep(DeBab, 1, FUN='*', y*1/DeBa^2))%*%(D%*%eBab_j)
    # second term
    eBab_jk <- sweep(B, 1, FUN='*', eBab_j)
    DeBab_jk <- D%*%eBab_jk
    second_term <- colSums(sweep(DeBab_jk, 1, FUN='*', y*(1/DeBa)))
    # third term
    third_term <- colSums(DeBab_jk)
    result[j,] <- first_term-second_term+third_term
  }
  return(result)
}
hess_a_negbin <- function(y, D, B, a, theta) {
  eBa <- exp(B%*%a)
  DeBa <- D%*%eBa
  #DeBa <- DeBa+1e-9
  eBab <- sweep(B, 1, FUN='*', eBa)
  DeBab <- D%*%eBab
  # fix j
  result <- matrix(NA, nrow=length(a), ncol=length(a))
  for (j in 1:length(a)) {
    eBab_j <- eBa*B[,j]
    eBab_jk <- sweep(B, 1, FUN='*', eBab_j)
    DeBab_jk <- D%*%eBab_jk
    first_term <- t(sweep(DeBab, 1, FUN='*', theta*1/DeBa^2))%*%(D%*%eBab_j)
    second_term <- colSums(sweep(DeBab_jk, 1, FUN='*', theta*(1/DeBa)))
    third_term <- t(sweep(DeBab, 1, FUN='*', y*1/(DeBa+theta)^2))%*%(D%*%eBab_j)
    fourth_term <- colSums(sweep(DeBab_jk, 1, FUN='*', y*(1/(DeBa+theta))))
    fifth_term <-  t(sweep(DeBab, 1, FUN='*', y*1/(DeBa)^2))%*%(D%*%eBab_j)
    sixth_term <- colSums(sweep(DeBab_jk, 1, FUN='*', y*(1/(DeBa))))
    result[j,] <- -first_term+second_term-third_term+fourth_term+fifth_term-sixth_term
  }
  return(result)
}

hess_a_theta_negbin <- function(y, D, B, a, theta) {
  eBa <- exp(B%*%a)
  DeBa <- D%*%eBa
  #DeBa <- DeBa+1e-9
  eBab <- sweep(B, 1, FUN='*', eBa)
  DeBab <- D%*%eBab
  n <- length(y)
  # fix j
  result <- matrix(NA, nrow=length(a)+1, ncol=length(a)+1)
  for (j in 1:length(a)) {
    eBab_j <- eBa*B[,j]
    eBab_jk <- sweep(B, 1, FUN='*', eBab_j)
    DeBab_jk <- D%*%eBab_jk
    first_term <- t(sweep(DeBab, 1, FUN='*', theta*1/(DeBa+theta)^2))%*%(D%*%eBab_j)
    second_term <- colSums(sweep(DeBab_jk, 1, FUN='*', theta*(1/(DeBa+theta))))
    third_term <- t(sweep(DeBab, 1, FUN='*', y*1/(DeBa)^2))%*%(D%*%eBab_j)
    fourth_term <- colSums(sweep(DeBab_jk, 1, FUN='*', y*(1/(DeBa))))
    fifth_term <-  t(sweep(DeBab, 1, FUN='*', y*1/(DeBa+theta)^2))%*%(D%*%eBab_j)
    sixth_term <- colSums(sweep(DeBab_jk, 1, FUN='*', y*(1/(DeBa+theta))))
    result[j,1:length(a)] <- -first_term+second_term-third_term+fourth_term+fifth_term-sixth_term
  }
  #dtheta^2
  dtheta2 <- sum(-trigamma(y+theta) + trigamma(theta) - (1/theta) +
                   2/(DeBa+theta) + (y+theta)*(DeBa+theta)^(-2))
  # dthetadalphaj
  dtheta_dalpha <- colSums(sweep(DeBab, 1, FUN='*', -(DeBa+theta)^(-1))) - 
    colSums(sweep(DeBab, 1, FUN='*', (y+theta)*(1/(DeBa+theta)^2)))
  # print(dtheta2)
  # print(dtheta_dalpha)
  result[length(a)+1,] <- c(dtheta_dalpha, dtheta2)
  result[,length(a)+1] <- c(dtheta_dalpha, dtheta2)
  return(result)
}

myquadprog_negbin <- function(a0, theta, y, D, B, conv = 0.001, maxiter=1000,
                              verbose=F) {
  lam_threshold = 1e-8; init_val <- -5; MIN_ITERATIONS <- 6
  beta_succ <- 1.1; beta_fail <- 0.5; gamma <- 0.1
  epsilon_2 <- 1e-5; delta <- 0.1 #trust region radius
  DF <- ncol(B)
  current_a <- a0
  current_theta <- theta
  loglik_memory <- seq(1,10)
  prev_loglik <- loglik_negbin(current_a,theta,y,D,B)
  diff <- 100
  iter <- 0
  while (iter < maxiter) {
    iter <- iter+1
    prev_a <- current_a
    prev_theta <- current_theta
    H_0 <- hess_a_theta_negbin(y,D,B,current_a,current_theta)
    grad_0 <- c(deriv_a_negbin(y,D,B,current_a,current_theta), 
                deriv_theta_negbin(y,D,B,current_a,current_theta))
    H <- psd(H_0)
    norm_factor <- norm(H,'2')
    H <- H/norm_factor
    grad <- grad_0/norm_factor
    epsilon <- 1e-7
    H <- H + epsilon*diag(length(grad))
    # grad <- grad_0
    A <- cbind(diag(length(grad)), -diag(length(grad)))
    bzero <- rep(-delta,2*dim(H)[2])
    solution <- quadprog::solve.QP(H,grad,A,bzero,meq=0)$solution
    # solution[1:DF] <- 0
    # solution[1+DF] <- 0
    # solution <- solution / 100
    sol_a <- solution[1:DF]
    sol_theta <- solution[DF+1]
    predicted_decrease <- as.numeric(-(0.5*t(solution)%*%H_0%*%solution - sum(grad_0*solution)))
    #print(solution)
    current_loglik <- loglik_negbin(current_a+sol_a,current_theta+sol_theta,y,D,B)
    true_decrease <- prev_loglik-current_loglik
    if (verbose) {
      print(paste('Previous log-likelihood:', prev_loglik))
      print(paste('Current log-likelihood:', current_loglik))
      print(paste('Predicted decrease:',predicted_decrease))
      print(paste('True decrease:',true_decrease))
    }
    myconst <- gamma*predicted_decrease
    if (verbose) {
      print(paste('Constant:', myconst))
      print(paste('True decrease >= constant?',true_decrease >= myconst))
    }
    diff <- prev_loglik-current_loglik
    if (true_decrease >= myconst) {
      # then we expand the trust region
      delta <- beta_succ*delta
      current_a <- current_a+sol_a
      current_theta <- current_theta+sol_theta
      prev_loglik <- current_loglik
    } else {
      # otherwise decrease the trust region
      delta <- min(1, beta_fail*delta)
    }
    loglik_memory <- c(loglik_memory, diff)
    if (length(loglik_memory)>10) {
      loglik_memory <- loglik_memory[-1]
    }
    if (verbose) {
      message(paste(diff))
      #message(paste(loglik_memory))
    }
    if (diff<conv) {
      if(abs(max(loglik_memory)-min(loglik_memory))<conv) {
        break
      }
    }
  }
  if (iter==maxiter) {
    message(paste('reached max iterations of',maxiter))
  }
  return(c(current_a,current_theta))
}
