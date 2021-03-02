###################################################################
# Apply to XD's package ####
source("/Volumes/GoogleDrive/My Drive/Work/HIV/straweib/R/icweib.R")
ci_SEF <- function(param_est, param_se, alpha = 0.05){
  c <- qnorm((1 - alpha/2))
  lb <- param_est - c*(param_se)
  ub <- param_est + c*(param_se)
  list(lb=lb, ub=ub)
}
xiangdongF <- function(n, a0, aj, beta, gamma, lambda){
  out <- simDataF(n, a0, aj, beta, gamma, lambda)
  L <- out$L
  R <- out$R
  eventtime <- out$t
  z <- out$z
  ## Fit weibull PH model with mother ARV only
  out <- icweib(L, R, data = out, covariates = ~z)
  beta_hat <- out$coef[1]
  beta_se <- out$coef[2]
  beta_ci <- ci_SEF(beta_hat, beta_se, alpha = 0.05)
  beta_logit <- beta_ci$lb < beta & beta < beta_ci$ub
  gamma_hat <- out$weib[3]
  gamma_se <- sqrt(out$cov[1,1])
  gamma_ci <- ci_SEF(gamma_hat, gamma_se, alpha = 0.05)
  gamma_logit <- gamma_ci$lb < gamma & gamma < gamma_ci$ub
  lambda_hat <- out$weib[4]
  lambda_se <- sqrt(out$cov[2,2])
  lambda_ci <- ci_SEF(lambda_hat, lambda_se, alpha = 0.05)
  lambda_logit <- lambda_ci$lb < lambda & lambda < lambda_ci$ub
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
gamma= 0.290028
lambda = 0.910090
out_0.9_pi0.1 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
tb_0.9_pi0.1 <- data.frame(beta_mean = mean(out_0.9_pi0.1[1,]),
                 beta_se = mean(out_0.9_pi0.1[2,]),
                 beta_coverage = mean(out_0.9_pi0.1[5,]),
                 gamma_mean = mean(out_0.9_pi0.1[6,]),
                 gamma_se = mean(out_0.9_pi0.1[7,]),
                 gamma_coverage = mean(out_0.9_pi0.1[10,]),
                 lambda_mean = mean(out_0.9_pi0.1[11,]),
                 lambda_se = mean(out_0.9_pi0.1[12,]),
                 lambda_coverage = mean(out_0.9_pi0.1[15,]))
write.csv(tb_0.9_pi0.1, 'Gu_tb_0.9_pi0.1.csv')
# pi=0.2#######################
a0 = -0.1
aj = -4
gamma=0.550054 
lambda=0.450044
out_0.9_pi0.2 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))

tb_0.9_pi0.2 <- data.frame(beta_mean = mean(out_0.9_pi0.2[1,]),
                           beta_se = mean(out_0.9_pi0.2[2,]),
                           beta_coverage = mean(out_0.9_pi0.2[5,]),
                           gamma_mean = mean(out_0.9_pi0.2[6,]),
                           gamma_se = mean(out_0.9_pi0.2[7,]),
                           gamma_coverage = mean(out_0.9_pi0.2[10,]),
                           lambda_mean = mean(out_0.9_pi0.2[11,]),
                           lambda_se = mean(out_0.9_pi0.2[12,]),
                           lambda_coverage = mean(out_0.9_pi0.2[15,]))
write.csv(tb_0.9_pi0.2, 'Gu_tb_0.9_pi0.2.csv')

# pi=0.4#######################
a0 = 1.6
aj = -4.4
gamma=3.000200e-02
lambda=7.100700e-01
out_0.9_pi0.4 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))

tb_0.9_pi0.4 <- data.frame(beta_mean = mean(out_0.9_pi0.4[1,]),
                           beta_se = mean(out_0.9_pi0.4[2,]),
                           beta_coverage = mean(out_0.9_pi0.4[5,]),
                           gamma_mean = mean(out_0.9_pi0.4[6,]),
                           gamma_se = mean(out_0.9_pi0.4[7,]),
                           gamma_coverage = mean(out_0.9_pi0.4[10,]),
                           lambda_mean = mean(out_0.9_pi0.4[11,]),
                           lambda_se = mean(out_0.9_pi0.4[12,]),
                           lambda_coverage = mean(out_0.9_pi0.4[15,]))
write.csv(tb_0.9_pi0.4, 'Gu_tb_0.9_pi0.4.csv')

################################### Failure Rate = 0.5 ##############################
failure_rate=0.5
# pi=0.1#######################
a0 = -0.9
aj = -4
gamma=0.690068
lambda=0.070006
out_0.5_pi0.1 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))

tb_0.5_pi0.1 <- data.frame(beta_mean = mean(out_0.5_pi0.1[1,]),
                           beta_se = mean(out_0.5_pi0.1[2,]),
                           beta_coverage = mean(out_0.5_pi0.1[5,]),
                           gamma_mean = mean(out_0.5_pi0.1[6,]),
                           gamma_se = mean(out_0.5_pi0.1[7,]),
                           gamma_coverage = mean(out_0.5_pi0.1[10,]),
                           lambda_mean = mean(out_0.5_pi0.1[11,]),
                           lambda_se = mean(out_0.5_pi0.1[12,]),
                           lambda_coverage = mean(out_0.5_pi0.1[15,]))
write.csv(tb_0.5_pi0.1, 'Gu_tb_0.5_pi0.1.csv')

# pi=0.2#######################
a0 = -0.1
aj = -4
gamma=0.570056
lambda=0.070006
out_0.5_pi0.2 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))
## beta
mean(out_0.5_pi0.2[1,])
mean(out_0.5_pi0.2[2,])
mean(out_0.5_pi0.2[5,])
## gamma
mean(out_0.5_pi0.2[6,])
mean(out_0.5_pi0.2[7,])
mean(out_0.5_pi0.2[10,])
## lambda
mean(out_0.5_pi0.2[11,])
mean(out_0.5_pi0.2[12,])
mean(out_0.5_pi0.2[15,])

tb_0.5_pi0.2 <- data.frame(beta_mean = mean(out_0.5_pi0.2[1,]),
                           beta_se = mean(out_0.5_pi0.2[2,]),
                           beta_coverage = mean(out_0.5_pi0.2[5,]),
                           gamma_mean = mean(out_0.5_pi0.2[6,]),
                           gamma_se = mean(out_0.5_pi0.2[7,]),
                           gamma_coverage = mean(out_0.5_pi0.2[10,]),
                           lambda_mean = mean(out_0.5_pi0.2[11,]),
                           lambda_se = mean(out_0.5_pi0.2[12,]),
                           lambda_coverage = mean(out_0.5_pi0.2[15,]))
write.csv(tb_0.5_pi0.2, 'Gu_tb_0.5_pi0.2.csv')

# pi=0.4#######################
a0 = 1.7
aj = -4.5
gamma = 2.500240e-01 
lambda = 5.000400e-02
out_0.5_pi0.4 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))

tb_0.5_pi0.4 <- data.frame(beta_mean = mean(out_0.5_pi0.4[1,]),
                           beta_se = mean(out_0.5_pi0.4[2,]),
                           beta_coverage = mean(out_0.5_pi0.4[5,]),
                           gamma_mean = mean(out_0.5_pi0.4[6,]),
                           gamma_se = mean(out_0.5_pi0.4[7,]),
                           gamma_coverage = mean(out_0.5_pi0.4[10,]),
                           lambda_mean = mean(out_0.5_pi0.4[11,]),
                           lambda_se = mean(out_0.5_pi0.4[12,]),
                           lambda_coverage = mean(out_0.5_pi0.4[15,]))

write.csv(tb_0.5_pi0.4, 'Gu_tb_0.5_pi0.4.csv')

################################### Failure Rate = 0.3 ##############################
failure_rate=0.3
# pi=0.1#######################
pi=0.1
gamma=0.530052
lambda=0.050004
a0 = -0.9
aj = -4
out_0.3_pi0.1 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))

tb_0.3_pi0.1 <- data.frame(beta_mean = mean(out_0.3_pi0.1[1,]),
                           beta_se = mean(out_0.3_pi0.1[2,]),
                           beta_coverage = mean(out_0.3_pi0.1[5,]),
                           gamma_mean = mean(out_0.3_pi0.1[6,]),
                           gamma_se = mean(out_0.3_pi0.1[7,]),
                           gamma_coverage = mean(out_0.3_pi0.1[10,]),
                           lambda_mean = mean(out_0.3_pi0.1[11,]),
                           lambda_se = mean(out_0.3_pi0.1[12,]),
                           lambda_coverage = mean(out_0.3_pi0.1[15,]))
write.csv(tb_0.3_pi0.1, 'Gu_tb_0.3_pi0.1.csv')

# pi=0.2#######################
pi=0.2
a0 = -4.2
aj = 4
gamma=3.300320e-01
lambda=5.000400e-02
out_0.3_pi0.2 <- replicate(nsim, xiangdongF(n, beta, a0, aj, gamma, lambda))

tb_0.3_pi0.2 <- data.frame(beta_mean = mean(out_0.3_pi0.2[1,]),
                           beta_se = mean(out_0.3_pi0.2[2,]),
                           beta_coverage = mean(out_0.3_pi0.2[5,]),
                           gamma_mean = mean(out_0.3_pi0.2[6,]),
                           gamma_se = mean(out_0.3_pi0.2[7,]),
                           gamma_coverage = mean(out_0.3_pi0.2[10,]),
                           lambda_mean = mean(out_0.3_pi0.2[11,]),
                           lambda_se = mean(out_0.3_pi0.2[12,]),
                           lambda_coverage = mean(out_0.3_pi0.2[15,]))
write.csv(tb_0.3_pi0.2, 'Gu_tb_0.3_pi0.2.csv')
