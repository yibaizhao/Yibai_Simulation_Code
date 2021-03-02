source('SimulationII_Hessian.R')
# Creat dataset
# Inverse CDF of weibull
## Function to generate pi
piF <- function(a0, aj, cov_aj){ # z is an indicator of treatment group
  if(is.null(cov_aj)){
    linear <- a0 + matrix(0, nrow = nrow(data))
  }else{
    linear <- a0 + cov_aj %*% aj
  }
  pi <- exp(linear) / (1 + exp(linear))
  return(pi)
}


eventtimeFc <- function(n, strata, cov_beta, pi, gamma, lambda, beta){ # gamma and lambda are parameters from Weibull(gamma, lambda)
  U <- runif(n) # Generate a uniformly distributed random variable U, in  order to find mixture T
  w <- runif(n) # Generate a uniformly distributed random variable w, representing F, where F=1-S
  eventtime <- c() 
  for (i in 1:n) {
    if(U[i]<=pi[i]){ # if U<=pi, 
      eventtime[i] <- 0 # P(T=0)=pi
    }else{ # if U>pi, P(T>0)=1-pi
      eventtime[i] <- (-log(1-w[i])/(lambda[which(strata[i,] == 1)]*exp(cov_beta[i,] %*% beta)))^(1/gamma[which(strata[i,] == 1)])
    }
  }
  return(eventtime)
}

simDataF <- function(n=n, a0, aj, beta, gamma, lambda){
  cov_beta <- as.matrix(sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))) # binary variable, like gender, cov_beta=0 as female, cov_beta=1 as male
  cov_aj <- strata <-  as.matrix(sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))) # Sample group indicator z with P(Z=0)=0.5
  if (ncol(cov_aj) == 1){
    strata <- cbind(1-cov_aj, cov_aj)
  }
  pi <- piF(a0, aj, cov_aj)
  k_init <-  sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5)) # define whether it's in mixture or not
  tt <- eventtimeFc(n, strata, cov_beta, pi, gamma, lambda, beta)
  tau <- c(1,4,10)
  k <- ifelse(tt>tau[1], 1, k_init) # k=1 if non-mixture, k=0 if mixture
  c <- ifelse(k==1 & tt == 0, 1,
              ifelse(k==1 & tt > 0, 0, NA))
    
  data <- data.frame(ID = rep(1:n, each=length(tau)), 
                     testtime = rep(tau, n),  
                     eventtime = rep(tt, each=length(tau)))
  data$result <- ifelse(data$eventtime >  data$testtime, 0, 1) # 0 if eventtime>testtime, 1 if not
  id <- unique(data$ID) # Find unique subject ID
  L0 <-
    sapply(id, function(x)
      tail(data[data$ID == x &
                  data$result == 0, ], 1)$testtime) # Find last negative to be left interval
  L <-
    sapply(1:length(id), function(x)
      ifelse(length(L0[[x]]) == 0, 0, L0[[x]])) # Set left interval to 0 if all are positive

  R0 <-
    sapply(id, function(x)
      head(data[data$ID == x &
                  data$result == 1, ], 1)$testtime) # Find first positive to be right interval
  R <-
    sapply(1:length(id), function(x)
      ifelse(length(R0[[x]]) == 0, Inf, R0[[x]]))# Set right interval to Inf if all are negative

  list(L=L, R=R, t = tt, strata=strata, cov_aj=cov_aj, cov_beta=cov_beta, c=c, k=k)
}

# Generate survival function under PH assumption when T~Weibull(gamma, lambda)
SurvF <- function(t, gamma, lambda, b, strata, cov_beta){
  gamma_all <- gamma[apply(strata, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(strata, 1, function(x) which(x == 1))]
  if(is.null(cov_beta)){
    exp(-lambda_all*t^gamma_all)
  }else{
    exp(-lambda_all*t^gamma_all*exp(cov_beta %*% b))
  }
}

# calculate ci
ciF <- function(param_est, param_var, alpha = 0.05){
  c <- qnorm((1 - alpha/2))
  lb <- param_est - c*sqrt(param_var)
  ub <- param_est + c*sqrt(param_var)
  list(lb=lb, ub=ub)
}


EiF <- function(beta, gamma, lambda, a0, aj, L, R, k, c, cov_aj, strata, cov_beta){
  pi <- piF(a0, aj, cov_aj)
  # in mixture model
  ei <- pi / (pi + (1-pi) * (1 - SurvF(R, gamma, lambda, beta, strata, cov_beta)))
  # Ei=1, in logistic model
  ei[k == 1 & c == 1] <- 1  
  # Ei=0, in survival model
  ei[k == 1 & c == 0] <- 0 
  ei
}
######################### EM Algorithm ######################
# Observed Loglikelihood ####
Loglik <- function(beta, gamma, lambda, a0, aj, L, R, cov_aj, strata, cov_beta, k, c){
  pi <- piF(a0, aj, cov_aj)
  SL <- SurvF(L, gamma, lambda, b=beta, strata, cov_beta)
  SR <- SurvF(R, gamma, lambda, b=beta, strata, cov_beta)
  
  loglik_i <- c * log(pi) + (1 - c) * log((1 - pi) * (SL - SR))
  loglik_i[k == 1 & c == 1] <- log(pi)[k == 1 & c == 1] # logit 
  loglik_i[k == 0] <- log(pi + (1 - pi) * (SL - SR))[k == 0] # mixture
  loglik <- sum(loglik_i)
  loglik ## times -1 to get maximum
  
}


# Expected Loglikelihood #####
ELoglik <- function(beta, gamma, lambda, a0, aj, L, R, cov_aj, strata, cov_beta, Ei){
  pi <- piF(a0, aj, cov_aj)
  SL <- SurvF(L, gamma, lambda, b=beta, strata, cov_beta)
  SR <- SurvF(R, gamma, lambda, b=beta, strata, cov_beta)

  loglik_i <- Ei * log(pi) + (1 - Ei) * log((1 - pi) * (SL - SR))
  loglik_i[SL == SR] <- (Ei * log(pi))[SL == SR]
  loglik <- sum(loglik_i)
  loglik ## times -1 to get maximum

}

parmsF <- function(nparms, L, R, cov_aj, strata, cov_beta, Ei){
  naj <- ncol(cov_aj)
  nbeta <- ncol(cov_beta)
  n_gamma_lambda <- ncol(strata)
  
  a0 <- nparms[1] 
  aj <- matrix(nparms[2 : (1+naj)], ncol = 1)
  beta <- matrix(nparms[(1+naj+1): (1+naj+nbeta)], ncol = 1)
  gamma <- exp(nparms[(1+naj+nbeta+1) : (1+naj+nbeta+n_gamma_lambda)])
  lambda <- exp(nparms[(1+naj+nbeta+n_gamma_lambda+1) : (1+naj+nbeta+n_gamma_lambda*2)])
  
  pi <- piF(a0, aj, cov_aj)
  SL <- SurvF(L, gamma, lambda, b=beta, strata, cov_beta)
  SR <- SurvF(R, gamma, lambda, b=beta, strata, cov_beta)
  
  loglik_i <- Ei * log(pi) + (1 - Ei) * log((1 - pi) * (SL - SR))
  loglik_i[SL == SR] <- (Ei * log(pi))[SL == SR]
  loglik <- sum(loglik_i)
  loglik
}

outF <- function(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0, total_num = 1000){
  num = 0 # counter for data points
  tb_ci <- c()
  tb_sd <- c()
  ests <- c()
  
  while (num < total_num) {
    ### Generate dataset
    out <- simDataF(n, a0_0, aj_0, beta0, gamma0, lambda0)
    L <- out$L
    R <- out$R
    cov_aj <- out$cov_aj
    strata <- out$strata
    cov_beta <- out$cov_beta
    k <- out$k
    c <- out$c
    
    naj <- ncol(cov_aj)
    nbeta <- ncol(cov_beta)
    n_gamma_lambda <- ncol(strata)
    
    ## Initialization
    ### Find initial values
    beta_init <- matrix(0.5, nrow = nbeta)
    aj_init <-  matrix(0.5, nrow = naj)
    a0_init <- 0.5
    log_gamma_init <- log_lambda_init <- rep(log(0.5), n_gamma_lambda)
    param_init <- c(a0_init, aj_init, beta_init, log_gamma_init, log_lambda_init)
    
    Ei_init <- c
    Ei_init[is.na(c)] <- sample(c(0, 1), sum(is.na(c)), replace = TRUE, prob = c(0.5, 0.5))
    est_nparms_init <- optim(par = param_init, parmsF, L=L, R=R, cov_aj=cov_aj, strata=strata, cov_beta=cov_beta, Ei=Ei_init,
                             control=list(fnscale=-1))$par # fnscale=-1 to find maximum
    
    a0 <- est_nparms_init[1] 
    aj <- matrix(est_nparms_init[2 : (1+naj)], ncol = 1)
    beta <- matrix(est_nparms_init[(1+naj+1): (1+naj+nbeta)], ncol = 1)
    log_gamma <- est_nparms_init[(1+naj+nbeta+1) : (1+naj+nbeta+n_gamma_lambda)]
    gamma <- exp(log_gamma)
    log_lambda <- est_nparms_init[(1+naj+nbeta+n_gamma_lambda+1) : (1+naj+nbeta+n_gamma_lambda*2)]
    lambda <- exp(log_lambda)
    est_nparms <- est_nparms_init
    
    ######################## End of Initialization ######################## 
    ## Maximization ####
    loglik <- 0
    loglik[2] <- ELoglik(beta, gamma, lambda, a0, aj, L, R, cov_aj, strata, cov_beta, Ei_init)

    i <- 2
    
    tol = 1e-6 # tolerance 1e-6
    while (abs(loglik[i] - loglik[i - 1]) >= tol) {
      # update Ei using theta^{t-1}
      Ei <- EiF(beta, gamma, lambda, a0, aj, L, R, k, c, cov_aj, strata, cov_beta)
      
      est_nparms <- optim(par = est_nparms, parmsF, L=L, R=R, cov_aj=cov_aj, strata=strata, cov_beta=cov_beta, Ei=Ei,
                               control=list(fnscale=-1))$par # fnscale=-1 to find maximum
      
      a0 <- est_nparms[1] 
      aj <- matrix(est_nparms[2 : (1+naj)], ncol = 1)
      beta <- matrix(est_nparms[(1+naj+1): (1+naj+nbeta)], ncol = 1)
      log_gamma <- est_nparms[(1+naj+nbeta+1) : (1+naj+nbeta+n_gamma_lambda)]
      gamma <- exp(log_gamma)
      log_lambda <- est_nparms[(1+naj+nbeta+n_gamma_lambda+1) : (1+naj+nbeta+n_gamma_lambda*2)]
      lambda <- exp(log_lambda)


      i <- i + 1
      loglik[i] <- ELoglik(beta, gamma, lambda, a0, aj, L, R, cov_aj, strata, cov_beta, Ei)
      # print(c(i, est_nparms, loglik[i]))
      print(paste0(i-2, "'s Iteration ..."))
    }
    print("Done!")
    est_nparms <- c(a0, aj, beta, gamma, lambda)
    maxLoglik <- Loglik(beta, gamma, lambda, a0, aj, L, R, cov_aj, strata, cov_beta, k, c)
    ests <- c(ests, i, est_nparms, maxLoglik)
    
    var_all <- try(diag(hessianM(a0, aj, beta, gamma, lambda, L, R, c, k, cov_aj, strata, cov_beta)), silent = TRUE)
    if ('try-error' %in% class(var_all)) next
    
    # calculate confidence interval and coverage probability
    ci_0.05 <- ciF(a0, var_all[1], alpha = 0.05)
    covLogit_a0 <- ci_0.05$lb < a0_0 & ci_0.05$ub > a0_0
    ci_0.05 <- ciF(aj, var_all[2 : (1+naj)], alpha = 0.05)
    covLogit_aj <- (ci_0.05$lb < aj_0 & ci_0.05$ub > aj_0)
    ci_0.05 <- ciF(beta, var_all[(1+naj+1): (1+naj+nbeta)], alpha = 0.05)
    covLogit_beta <- (ci_0.05$lb < beta0 & ci_0.05$ub > beta0)
    ci_0.05 <- ciF(gamma, var_all[(1+naj+nbeta+1) : (1+naj+nbeta+n_gamma_lambda)], alpha = 0.05)
    covLogit_gamma <- (ci_0.05$lb < gamma0 & ci_0.05$ub > gamma0)
    ci_0.05 <- ciF(lambda, var_all[(1+naj+nbeta+n_gamma_lambda+1) : (1+naj+nbeta+n_gamma_lambda*2)], alpha = 0.05)
    covLogit_lambda <- (ci_0.05$lb < lambda0 & ci_0.05$ub > lambda0)
    
    tb_ci <- c(tb_ci, covLogit_a0, covLogit_aj, covLogit_beta, covLogit_gamma, covLogit_lambda)
    tb_sd <- c(tb_sd, sqrt(var_all))
    
    num = num + 1
    print(paste0("# of Datapoints: ", num))
          
  }
  
  nparams <- (1+naj+nbeta+n_gamma_lambda*2)
  out_ests0 <- matrix(ests, byrow = TRUE, ncol = nparams+2)[,-c(1, nparams+2)]

  out_ci <- colMeans(matrix(tb_ci, byrow = TRUE, ncol = nparams))
  out_ests <- colMeans(out_ests0)
  out_sd <- colMeans(matrix(tb_sd, byrow = TRUE, ncol = nparams))
  
  out0 <- rbind(out_ests, out_sd, out_ci)
  out <- as.data.frame(t(out0))
  rownames(out) <- c('a0', 'aj', 'beta', paste0('gamma', 1:n_gamma_lambda), paste0('lambda', 1:n_gamma_lambda))
  colnames(out) <- c('Estimates', 'Std', 'CovProb')
  
  write.csv(out, file = paste0("Results/out_f", failure_rate, "_pi", pi0, '_', Sys.Date(), ".csv"))
  write.csv(out_ests0, file = paste0('Results/', failure_rate, "_pi", pi0, '_ests_', Sys.Date(), ".csv"))
  
}
