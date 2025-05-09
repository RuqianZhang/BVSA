setwd("./")
source("BVSA.R")
require(foreach)
require(doParallel)
require(MASS)

Coef_gen <- function(p, theta0, alpha0, s, rho){
  beta0 <- rep(0,p+1); beta0[1:(s+1)] <- theta0
  gamma0 <- rep(0,p+1); gamma0[1:(s+1)] <- theta0
  alpha0 <- alpha0
  
  return(list(beta0=beta0, gamma0=gamma0, alpha0=alpha0))
}

numCores <- max(detectCores()-3, 1)
registerDoParallel(numCores)

theta_list <- list(c(1, -1.5, 2, -2.5, 3))
alpha_list <- list(c(40, 0))
n_list <- seq(100, 1000, 100)


setwd("./simu/res/est")
for (p in c(50)){
  for (n in n_list){
    for (rho in c(0)){
      for (s in c(4)){
        for (theta in theta_list){
          for (alpha in alpha_list){
            cat("p:", p, "n:", n, "rho:", rho, "s:", s, "theta0:", theta, "alpha0:", alpha, '\n')
            
            # Generate parameters
            params <- Coef_gen(p, theta, alpha, s, rho)
            beta0 <- params$beta0; gamma0 <- params$gamma0; alpha0 = params$alpha0
            
            ntrial <- 100
            measures <- c("TP_b", "TP_s_b", "FP_b", "exact_b", "included_b", "exact_s_b", "l2_error_b", 
                          "TP_g", "TP_s_g", "FP_g", "exact_g", "included_g", "exact_s_g", "l2_error_g",
                          "l2_error_a", "subgroup_pred_error", "time", "seed")
            record_perf <- matrix(0, ncol = length(measures), nrow = ntrial+2)
            colnames(record_perf) <- measures
            
            record_gamma_estimate <- matrix(0, ncol = p+1, nrow = ntrial+2)
            record_beta_estimate <- matrix(0, ncol = p+1, nrow = ntrial+2)
            record_alpha_estimate <- matrix(0, ncol = 2, nrow = ntrial+2)
            colnames(record_alpha_estimate) <- c("alpha_1", "alpha_2")
            
            loops <- ceiling(ntrial/numCores)
            for (i in 1:loops){
              runs <- ifelse(i < loops, numCores, ntrial-(loops-1)*numCores)
              Loop_Res <- foreach(run=1:runs, .combine=rbind) %dopar% {
                # Generate sample
                run_seed <- 28+((i-1)*numCores+run)*2
                set.seed(run_seed)
                cov1 <- (1-rho)*diag(s) + array(rho, c(s, s))
                cov3 <- (1-rho)*diag(p-s) + array(rho, c(p-s, p-s))
                cov2 <- array(rho, c(s, p-s))
                cov_mat <- rbind(cbind(cov1, cov2), cbind(t(cov2), cov3))
                X <- cbind(rep(1,n), mvrnorm(n, rep(0,p), cov_mat))
                Z <- cbind(rep(1,n), mvrnorm(n, rep(0,p), cov_mat))
                X[,-1] <- scale(X[,-1])*sqrt(n/(n-1))
                Z[,-1] <- scale(Z[,-1])*sqrt(n/(n-1))
                X_types <- rep("s", p)
                
                trt <- rbinom(n,1,0.5)
                p_delta0 <- 1/(1+exp(-X%*%gamma0))
                delta0 <- rbinom(n, 1, p_delta0)
                sigmasq0 <- 1
                eps <- rnorm(n, 0, sigmasq0)
                Y <- Z%*%beta0+delta0*trt*alpha0[1]+(1-delta0)*trt*alpha0[2]+eps
                
                # Fitting
                tau_gamma0=5/n; q_gamma=min(0.2, 10/p); tau_beta0=5/n; q_beta=min(0.2, 10/p)
                tau_gamma1=max(sqrt(p^2/(100*n)),1)
                tau_beta1=max(sqrt(p^2/(100*n)),1)
                
                start <- Sys.time()
                result_BVSA <- BVSA_C(X, Y, trt, Z, iter=20000, burn=5000,
                                               tau_gamma0=tau_gamma0, q_gamma=q_gamma, tau_gamma1=tau_gamma1,
                                               tau_beta0=tau_beta0, q_beta=q_beta, tau_beta1=tau_beta1)
                end <- Sys.time()
                diff_time <- difftime(end, start, units="secs")[[1]]
                Igamma <- result_BVSA[[1]]; p_gamma <- result_BVSA[[2]]
                Ibeta <- result_BVSA[[3]]; p_beta <- result_BVSA[[4]]
                Igamma_fix <- fixed_size(p_gamma, s)
                Ibeta_fix <- fixed_size(p_beta, s)
                gamma <- result_BVSA[[6]]; beta <- result_BVSA[[7]]
                alpha <- result_BVSA[[5]]
                sub_prob <- 1/(1+exp(-as.matrix(X[,which(Igamma==1)])%*%as.matrix(gamma[which(Igamma==1)])))
                subgroup_error <- sum(abs((sub_prob>=0.5)-delta0))/n
                
                res <- rep(0,length(measures))
                res[c(1,3,4,5)] <- perf(cbind(Ibeta, beta0))
                res[c(2,6)] <- perf(cbind(Ibeta_fix, beta0))[c(1,3)]
                res[7] <- sqrt(sum((beta-beta0)^2))
                res[c(8,10,11,12)] <- perf(cbind(Igamma, gamma0))
                res[c(9,13)] <- perf(cbind(Igamma_fix, gamma0))[c(1,3)]
                res[14] <- sqrt(sum((gamma-gamma0)^2))
                res[15] <- sqrt(sum((alpha-alpha0)^2))
                res[16] <- subgroup_error
                res[17] <- diff_time
                res[18] <- run_seed
                
                list(res, gamma, beta, alpha)
              }
              stopImplicitCluster()
              for (run in 1:runs){
                record_perf[(i-1)*numCores+run,] <- Loop_Res[run,][[1]]
                record_gamma_estimate[(i-1)*numCores+run,] <- Loop_Res[run,][[2]]
                record_beta_estimate[(i-1)*numCores+run,] <- Loop_Res[run,][[3]]
                record_alpha_estimate[(i-1)*numCores+run,] <- Loop_Res[run,][[4]]
              }
              
              write.csv(record_perf, paste("p", p, " n", n, " rho", rho, " Skinny perf.csv", sep=''), row.names = F)
              write.csv(record_gamma_estimate, paste("p", p, " n", n, " rho", rho, " Skinny gamma estimate.csv", sep=''), row.names = F)
              write.csv(record_beta_estimate, paste("p", p, " n", n, " rho", rho, " Skinny beta estimate.csv", sep=''), row.names = F)
              write.csv(record_alpha_estimate, paste("p", p, " n", n, " rho", rho, " Skinny alpha estimate.csv", sep=''), row.names = F)
              
              cat("Loop ", i, "done.", '\n', sep = '')
            }
            record_perf[ntrial+1,] <- apply(record_perf[1:ntrial,], 2, mean)
            record_gamma_estimate[ntrial+1,] <- apply(record_gamma_estimate[1:ntrial,], 2, mean)
            record_beta_estimate[ntrial+1,] <- apply(record_beta_estimate[1:ntrial,], 2, mean)
            record_alpha_estimate[ntrial+1,] <- apply(record_alpha_estimate[1:ntrial,], 2, mean)
            
            record_perf[ntrial+2,] <- apply(record_perf[1:ntrial,], 2, sd)
            record_gamma_estimate[ntrial+2,] <- apply(record_gamma_estimate[1:ntrial,], 2, sd)
            record_beta_estimate[ntrial+2,] <- apply(record_beta_estimate[1:ntrial,], 2, sd)
            record_alpha_estimate[ntrial+2,] <- apply(record_alpha_estimate[1:ntrial,], 2, sd)
            
            write.csv(record_perf, paste("p", p, " n", n, " rho", rho, " Skinny perf.csv", sep=''), row.names = F)
            write.csv(record_gamma_estimate, paste("p", p, " n", n, " rho", rho, " Skinny gamma estimate.csv", sep=''), row.names = F)
            write.csv(record_beta_estimate, paste("p", p, " n", n, " rho", rho, " Skinny beta estimate.csv", sep=''), row.names = F)
            write.csv(record_alpha_estimate, paste("p", p, " n", n, " rho", rho, " Skinny alpha estimate.csv", sep=''), row.names = F)
          }
        }
      }
    }
  }
}

stopImplicitCluster()















