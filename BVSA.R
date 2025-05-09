library(Rcpp)
sourceCpp("BVSA_main.cpp")

fixed_size <- function(theta, sparsity){
  max_theta <- order(theta[-1], decreasing=T)[1:sparsity]+1
  Itheta <- c(1, rep(0, length(theta)-1)); Itheta[max_theta] <- 1
  return(Itheta)
}

perf <- function(theta_pair){
  theta <- theta_pair[,1]
  theta0 <- theta_pair[,2]
  Itheta0 <- as.numeric(theta0[-1] != 0)
  Itheta <- as.numeric(theta[-1] != 0)
  TP <- sum(Itheta[which(Itheta0 == 1)] == 1)
  FP <- sum(Itheta[which(Itheta0 == 0)] == 1)
  exact <- as.numeric(all(Itheta == Itheta0))
  included <- as.numeric(TP == sum(Itheta0))
  return(c(TP, FP, exact, included))
}
