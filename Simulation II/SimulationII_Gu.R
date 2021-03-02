###################################################################
# Apply to XD's package ####
source("/Users/yibai/Google Drive/Work/HIV/straweib/R/icweib.R")
ci_SEF <- function(param_est, param_se, alpha = 0.05){
  c <- qnorm((1 - alpha/2))
  lb <- param_est - c*(param_se)
  ub <- param_est + c*(param_se)
  list(lb=lb, ub=ub)
}
xiangdongF <- function(n, beta, a0, aj, gamma, lambda){
  out <- simDataF(n, a0, aj, beta, gamma, lambda)
  L <- out$L
  R <- out$R
  eventtime <- out$t
  z <- out$z
  cov <- out$cov
  ## Fit weibull PH model with mother ARV only
  out <- icweib(L, R, data = out, strata = z, covariates = ~cov)
  # out <- icweib(L, R, data = out, covariates = ~z)
  beta_hat <- out$coef[1]
  beta_se <- out$coef[2]
  beta_ci <- ci_SEF(beta_hat, beta_se, alpha = 0.05)
  beta_logit <- beta_ci$lb < beta & beta < beta_ci$ub
  gamma_hat <- out$weib[,3]
  gamma_se <- sqrt(c(out$cov[1,1], out$cov[2,2]))
  gamma_ci <- ci_SEF(gamma_hat, gamma_se, alpha = 0.05)
  gamma_logit <- gamma_ci$lb < gamma & gamma < gamma_ci$ub
  lambda_hat <- out$weib[,4]
  lambda_se <- sqrt(c(out$cov[3,3], out$cov[4,4]))
  lambda_ci <- ci_SEF(lambda_hat, lambda_se, alpha = 0.05)
  lambda_logit <- lambda_ci$lb < lambda & lambda < lambda_ci$ub
  # beta_prob <- mean(beta_logit)
  c(as.numeric(beta_hat), as.numeric(beta_se), as.numeric(beta_ci$lb), as.numeric(beta_ci$ub), beta_logit,
    as.numeric(gamma_hat), as.numeric(gamma_se), as.numeric(gamma_ci$lb), as.numeric(gamma_ci$ub), gamma_logit,
    as.numeric(lambda_hat), as.numeric(lambda_se), as.numeric(lambda_ci$lb), as.numeric(lambda_ci$ub), lambda_logit)
}
nsim=1000 # number of iteration
set.seed(100)
beta = 0.7
n = 2500 # number of objects
################################### Failure Rate = 0.9 ##############################
failure_rate = 0.9
# pi=0.1#######################
a0 = -0.9
aj = -4
gamma = c(0.1, 0.2)
lambda = c(2, 1.2)
out_0.9_pi0.1 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.9_pi0.1_rowMean <- rowMeans(out_0.9_pi0.1)
tb_0.9_pi0.1 <- as.data.frame(
  as.matrix(tb_0.9_pi0.1_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.9_pi0.1, 'Gu_tb_0.9_pi0.1.csv')

# pi=0.2#######################
a0 = -0.6
aj = -3
gamma = c(0.1, 0.2)
lambda = c(2.3, 1.1)

out_0.9_pi0.2 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.9_pi0.2_rowMean <- rowMeans(out_0.9_pi0.2)
tb_0.9_pi0.2 <- as.data.frame(
  as.matrix(tb_0.9_pi0.2_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.9_pi0.2, 'Gu_tb_0.9_pi0.2.csv')

# pi=0.4#######################
a0 = 2.2
aj = -7
gamma = c(0.1, 0.3)
lambda = c(1.2, 0.9)

out_0.9_pi0.4 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.9_pi0.4_rowMean <- rowMeans(out_0.9_pi0.4)
tb_0.9_pi0.4 <- as.data.frame(
  as.matrix(tb_0.9_pi0.4_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.9_pi0.4, 'Gu_tb_0.9_pi0.4.csv')

################################### Failure Rate = 0.5 ##############################
failure_rate=0.5
# pi=0.1#######################
a0 = -0.9
aj = -4
gamma = c(0.1, 0.2)
lambda = c(1.0, 0.1)
out_0.5_pi0.1 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.5_pi0.1_rowMean <- rowMeans(out_0.5_pi0.1)
tb_0.5_pi0.1 <- as.data.frame(
  as.matrix(tb_0.5_pi0.1_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.5_pi0.1, 'Gu_tb_0.5_pi0.1.csv')

# pi=0.2#######################
a0 = -0.1
aj = -4
gamma = c(0.1, 0.2)
lambda = c(0.8, 0.1)
out_0.5_pi0.2 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.5_pi0.2_rowMean <- rowMeans(out_0.5_pi0.2)
tb_0.5_pi0.2 <- as.data.frame(
  as.matrix(tb_0.5_pi0.2_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.5_pi0.2, 'Gu_tb_0.5_pi0.2.csv')

# pi=0.4#######################
a0 = 1.5
aj = -4
gamma = c(0.3, 0.1)
lambda = c(0.2, 0.1)
out_0.5_pi0.4 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.5_pi0.4_rowMean <- rowMeans(out_0.5_pi0.4)
tb_0.5_pi0.4 <- as.data.frame(
  as.matrix(tb_0.5_pi0.4_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.5_pi0.4, 'Gu_tb_0.5_pi0.4.csv')


################################### Failure Rate = 0.3 ##############################
failure_rate=0.3
# pi=0.1#######################
a0 = -0.9
aj = -4
gamma = c(0.1, 0.2)
lambda = c(0.1, 0.1)
out_0.3_pi0.1 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.3_pi0.1_rowMean <- rowMeans(out_0.3_pi0.1)
tb_0.3_pi0.1 <- as.data.frame(
  as.matrix(tb_0.3_pi0.1_rowMean[c(1, 2, 5, 6, 8, 14, 7, 9, 15, 16, 18, 24, 17, 19, 25)]),
  row.names = c('beta_mean', 'beta_se', 'beta_coverage', 
                'gamma0_mean', 'gamma0_se', 'gamma0_coverage', 
                'gammaj_mean', 'gammaj_se', 'gammaj_coverage',
                'lambda0_mean', 'lambda0_se', 'lambda0_coverage',
                'lambdaj_mean', 'lambdaj_se', 'lambdaj_coverage')
)
write.csv(tb_0.3_pi0.1, 'Gu_tb_0.3_pi0.1.csv')