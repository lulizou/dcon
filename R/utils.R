
# Get gaussian basis functions over an interval (vector)
constructB <- function(interval, df = 10) {
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