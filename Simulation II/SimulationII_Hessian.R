######################################## Hessian matrix ######################################## 
exp_cov_beta <- function(beta, cov){
  if(any(is.na(beta))){
    1
  }else{
    exp(cov %*% beta)
  }
}
## logistic model ####
pi_dev_a0 <- function(pi){pi - pi^2}

pi_dev_aj <- function(pi, z_aj){z_aj * (pi - pi^2)}

pi_dev_a0_a0 <- function(pi){(1 - 2*pi) * pi_dev_a0(pi)}

pi_dev_aj_aj <- function(pi, z_aj, z_aj2){z_aj * (1 - 2*pi) * pi_dev_aj(pi, z_aj2)}

pi_dev_a0_aj <- function(pi, z_aj){(1 - 2*pi) * pi_dev_aj(pi, z_aj)}

## survival model ####
S_dev_lambda <- function(z, cov, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  SurvF(t, gamma, lambda, beta, z, cov) * (-t^gamma_all * exp_cov_beta(beta, cov))
}

S_dev_gamma <- function(z, cov, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  SurvF(t, gamma, lambda, beta, z, cov) * (-lambda_all*exp_cov_beta(beta, cov)) * t^gamma_all * log(t)
}

S_dev_beta <- function(z, cov, cov_j, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  SurvF(t, gamma, lambda, beta, z, cov) * (-lambda_all*t^gamma_all) * exp_cov_beta(beta, cov) * cov_j
}

S_dev_lambda_lambda <- function(z, cov, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  -t^gamma_all * exp_cov_beta(beta, cov) * S_dev_lambda(z, cov, t, beta, gamma, lambda)
}

S_dev_gamma_gamma <- function(z, cov, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  
  -lambda_all * exp_cov_beta(beta, cov) * log(t) * t^gamma_all *
    (S_dev_gamma(z, cov, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z, cov) * log(t))
}

S_dev_lambda_gamma <- function(z, cov, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  
  -exp_cov_beta(beta, cov) * t^gamma_all * (S_dev_gamma(z, cov, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z, cov)*log(t))
}

S_dev_lambda_beta <- function(z, cov, cov_j, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  
  -t^gamma_all * exp_cov_beta(beta, cov) * (S_dev_beta(z, cov, cov_j, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z, cov)*cov_j)
}

S_dev_gamma_beta <- function(z, cov, cov_j, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  
  -lambda_all * t^gamma_all * log(t) * exp_cov_beta(beta, cov) * (S_dev_beta(z, cov, cov_j, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z, cov)*cov_j)
}

S_dev_beta_beta <- function(z, cov, cov_j, cov_j2, t, beta, gamma, lambda){
  gamma_all <- gamma[apply(z, 1, function(x) which(x == 1))]
  lambda_all <- lambda[apply(z, 1, function(x) which(x == 1))]
  
  -lambda_all * t^gamma_all * cov_j * exp_cov_beta(beta, cov) * (S_dev_beta(z, cov, cov_j2, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z, cov)*cov_j2)
}

## derivative of Q(theta) ####
Q_dev_a0 <- function(pi, Ei){
  Ei * (1/pi) * pi_dev_a0(pi) - (1 - Ei) * (1 / (1 - pi)) * pi_dev_a0(pi)
}

Q_dev_aj <- function(pi, Ei, z_aj){
  Ei * (1/pi) * pi_dev_aj(pi, z_aj) - (1 - Ei) * (1 / (1 - pi)) * pi_dev_aj(pi, z_aj)
}

Q_dev_a0_a0 <- function(pi, Ei){
  Ei * (- (1 / pi^2) * pi_dev_a0(pi)^2 + (1 / pi) * pi_dev_a0_a0(pi)) - 
    (1 - Ei) * ((1 / (1 - pi)^2) * pi_dev_a0(pi)^2 + (1 / (1 - pi)) * pi_dev_a0_a0(pi))
}

Q_dev_aj_aj <- function(pi, Ei, z_aj, z_aj2){
  Ei * (- (1/pi^2) * pi_dev_aj(pi, z_aj) * pi_dev_aj(pi, z_aj2) + 1/pi * pi_dev_aj_aj(pi, z_aj, z_aj2)) -
    (1 - Ei) * ((1 / (1 - pi)^2) * pi_dev_aj(pi, z_aj) * pi_dev_aj(pi, z_aj2) + (1 / (1 - pi)) * pi_dev_aj_aj(pi, z_aj, z_aj2))
}

Q_dev_a0_aj <- function(pi, Ei, z_aj){
  Ei * (- (1 / pi^2) * pi_dev_aj(pi, z_aj) * pi_dev_a0(pi) + (1/ pi) * pi_dev_a0_aj(pi, z_aj)) -
    (1 - Ei) * ((1 / (1 - pi)^2) * pi_dev_aj(pi, z_aj) * pi_dev_a0(pi) + (1 / (1-pi)) * pi_dev_a0_aj(pi, z_aj))
}

Q_dev_beta <- function(Ei, z, cov, cov_j, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  q_dev_b <- (1 - Ei) * (1 / (SL - SR)) * (SL_dev_beta - SR_dev_beta)
  q_dev_b[Ei == 1] <- 0
  q_dev_b
}

Q_dev_lambda <- function(Ei, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  q_dev_lambda <- (1 - Ei) * (1 / (SL - SR)) * (SL_dev_lambda - SR_dev_lambda)
  q_dev_lambda[Ei == 1] <- 0
  q_dev_lambda
  
}

Q_dev_gamma <- function(Ei, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  q_dev_gamma <- (1 - Ei) * (1 / (SL - SR)) * (SL_dev_gamma - SR_dev_gamma)
  q_dev_gamma[Ei == 1] <- 0
  q_dev_gamma
}

Q_dev_beta_beta <- function(Ei, z, cov, cov_j, cov_j2, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_beta2 <- S_dev_beta(z, cov, cov_j2, L, beta, gamma, lambda)
  SL_dev_beta2[L == 0] <- 0
  SR_dev_beta2 <- S_dev_beta(z, cov, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  
  SL_dev_beta_beta <- S_dev_beta_beta(z, cov, cov_j, cov_j2, L, beta, gamma, lambda)
  SL_dev_beta_beta[L == 0] <- 0
  SR_dev_beta_beta <- S_dev_beta_beta(z, cov, cov_j, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta_beta[R == Inf] <- 0
  
  q_dev_b_b <- (1 - Ei) * (- (1 / (SL - SR)^2) * (SL_dev_beta - SR_dev_beta) * (SL_dev_beta2 - SR_dev_beta2) +
                             1 / (SL - SR) * (SL_dev_beta_beta - SR_dev_beta_beta)
  )
  q_dev_b_b[Ei == 1] <- 0
  q_dev_b_b
}

Q_dev_lambda_lambda <- function(Ei, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_lambda <- S_dev_lambda_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda_lambda[L == 0] <- 0
  SR_dev_lambda_lambda <- S_dev_lambda_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda_lambda[R == Inf] <- 0
  
  
  q_dev_lambda_lambda <- (1 - Ei) * (- (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda)^2 +
                                       1 / (SL - SR) * (SL_dev_lambda_lambda - SR_dev_lambda_lambda))
  q_dev_lambda_lambda[Ei == 1] <- 0
  q_dev_lambda_lambda
  
}

Q_dev_gamma_gamma <- function(Ei, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_gamma <- S_dev_gamma_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma_gamma[L == 0] <- 0
  SR_dev_gamma_gamma <- S_dev_gamma_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma_gamma[R == Inf] <- 0
  
  q_dev_gamma_gamma <- (1 - Ei) * (- (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma)^2 +
                                     1 / (SL - SR) * (SL_dev_gamma_gamma - SR_dev_gamma_gamma))
  q_dev_gamma_gamma[Ei == 1] <- 0
  q_dev_gamma_gamma
  
}

Q_dev_lambda_beta <- function(Ei, z, cov, cov_j, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_lambda_beta <- S_dev_lambda_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_lambda_beta[L == 0] <- 0
  SR_dev_lambda_beta <- S_dev_lambda_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_lambda_beta[R == Inf] <- 0
  
  q_dev_lambda_beta <- (1 - Ei) * (
    - (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_lambda_beta - SR_dev_lambda_beta)
  ) 
  q_dev_lambda_beta[Ei == 1] <- 0
  q_dev_lambda_beta
  
}

Q_dev_gamma_beta <- function(Ei, z, cov, cov_j, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_beta <- S_dev_gamma_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_gamma_beta[L == 0] <- 0
  SR_dev_gamma_beta <- S_dev_gamma_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_gamma_beta[R == Inf] <- 0
  
  q_dev_gamma_beta <- (1 - Ei) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_gamma_beta - SR_dev_gamma_beta)
  ) 
  q_dev_gamma_beta[Ei == 1] <- 0
  q_dev_gamma_beta
  
}

Q_dev_lambda_gamma <- function(Ei, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_gamma <- S_dev_lambda_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda_gamma[L == 0] <- 0
  SR_dev_lambda_gamma <- S_dev_lambda_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda_gamma[R == Inf] <- 0
  
  
  q_dev_lambda_gamma <- (1 - Ei) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_lambda - SR_dev_lambda) +
      1 / (SL - SR) * (SL_dev_lambda_gamma - SR_dev_lambda_gamma)
  ) 
  q_dev_lambda_gamma[Ei == 1] <- 0
  q_dev_lambda_gamma
  
}

## Observed loglikelihood for K0 group
k0_dev_a0_a0 <- function(mx, pi, beta, gamma, lambda, R, z, cov){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  (1 / mx) * SR * (pi_dev_a0_a0(pi) - (1 / mx) * SR * pi_dev_a0(pi)^2)
}

k0_dev_aj_aj <- function(mx, pi, beta, gamma, lambda, R, z, cov, z_aj, z_aj2){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  SR * (1 / mx) * (pi_dev_aj_aj(pi, z_aj, z_aj2) - (1 / mx) * SR * pi_dev_aj(pi, z_aj) * pi_dev_aj(pi, z_aj2))
}

k0_dev_a0_aj <- function(mx, pi, beta, gamma, lambda, R, z, cov, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  SR * (1/mx) * (pi_dev_a0_aj(pi, z_aj) - (1/mx) * SR * pi_dev_a0(pi) * pi_dev_aj(pi, z_aj))
}

k0_dev_beta_beta <- function(mx, pi, beta, gamma, lambda, R, z, cov, cov_j, cov_j2){
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  SR_dev_beta2 <- S_dev_beta(z, cov, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  SR_dev_beta_beta <- S_dev_beta_beta(z, cov, cov_j, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta_beta[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx)*(1-pi) * SR_dev_beta * SR_dev_beta2 + SR_dev_beta_beta)
}

k0_dev_gamma_gamma <- function(mx, pi, beta, gamma, lambda, R, z, cov){
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  SR_dev_gamma_gamma <- S_dev_gamma_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma_gamma[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_gamma^2 + SR_dev_gamma_gamma)
}

k0_dev_lambda_lambda <- function(mx, pi, beta, gamma, lambda, R, z, cov){
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  SR_dev_lambda_lambda <- S_dev_lambda_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda_lambda[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_lambda^2 + SR_dev_lambda_lambda)
}

k0_dev_gamma_beta <- function(mx, pi, beta, gamma, lambda, R, z, cov, cov_j){
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  SR_dev_gamma_beta <- S_dev_gamma_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_gamma_beta[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_gamma * SR_dev_beta + SR_dev_gamma_beta)
}

k0_dev_lambda_beta <- function(mx, pi, beta, gamma, lambda, R, z, cov, cov_j){
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  SR_dev_lambda_beta <- S_dev_lambda_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_lambda_beta[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_lambda * SR_dev_beta + SR_dev_lambda_beta)
}

k0_dev_lambda_gamma <- function(mx, pi, beta, gamma, lambda, R, z, cov){
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  SR_dev_lambda_gamma <- S_dev_lambda_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda_gamma[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_lambda * SR_dev_gamma + SR_dev_lambda_gamma)
}

k0_dev_a0_beta <- function(mx, pi, beta, gamma, lambda, R, z, cov, cov_j){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  pi_dev_a0(pi) * (1 / mx) * SR_dev_beta * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_a0_gamma <- function(mx, pi, beta, gamma, lambda, R, z, cov){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  pi_dev_a0(pi) * (1 / mx) * SR_dev_gamma * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_a0_lambda <- function(mx, pi, beta, gamma, lambda, R, z, cov){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  pi_dev_a0(pi) * (1 / mx) * SR_dev_lambda * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_aj_beta <- function(mx, pi, beta, gamma, lambda, R, z, cov, cov_j2, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SR_dev_beta2 <- S_dev_beta(z, cov, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  
  pi_dev_aj(pi, z_aj) * (1 / mx) * SR_dev_beta2 * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_aj_gamma <- function(mx, pi, beta, gamma, lambda, R, z, cov, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  pi_dev_aj(pi, z_aj) * (1 / mx) * SR_dev_gamma * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_aj_lambda <- function(mx, pi, beta, gamma, lambda, R, z, cov, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  pi_dev_aj(pi, z_aj) * (1 / mx) * SR_dev_lambda * (1 + (1 / mx) * (1 - pi) * SR)
}

k1_dev_a0_a0 <- function(pi, c){
  c * (- (1 / pi^2) * pi_dev_a0(pi)^2 + (1 / pi) * pi_dev_a0_a0(pi)) - 
    (1 - c) * ((1 / (1 - pi)^2) * pi_dev_a0(pi)^2 + (1 / (1 - pi)) * pi_dev_a0_a0(pi))
}

k1_dev_aj_aj <- function(pi, c, z_aj, z_aj2){
  c * (- (1/pi^2) * pi_dev_aj(pi, z_aj) * pi_dev_aj(pi, z_aj2) + 1/pi * pi_dev_aj_aj(pi, z_aj, z_aj2)) -
    (1 - c) * ((1 / (1 - pi)^2) * pi_dev_aj(pi, z_aj) * pi_dev_aj(pi, z_aj2) + (1 / (1 - pi)) * pi_dev_aj_aj(pi, z_aj, z_aj2))
}

k1_dev_a0_aj <- function(pi, c, z_aj){
  c * (- (1 / pi^2) * pi_dev_aj(pi, z_aj) * pi_dev_a0(pi) + (1/ pi) * pi_dev_a0_aj(pi, z_aj)) -
    (1 - c) * ((1 / (1 - pi)^2) * pi_dev_aj(pi, z_aj) * pi_dev_a0(pi) + (1 / (1-pi)) * pi_dev_a0_aj(pi, z_aj))
}

k1_dev_beta_beta <- function(c, z, cov, cov_j, cov_j2, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_beta2 <- S_dev_beta(z, cov, cov_j2, L, beta, gamma, lambda)
  SL_dev_beta2[L == 0] <- 0
  SR_dev_beta2 <- S_dev_beta(z, cov, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  
  SL_dev_beta_beta <- S_dev_beta_beta(z, cov, cov_j, cov_j2, L, beta, gamma, lambda)
  SL_dev_beta_beta[L == 0] <- 0
  SR_dev_beta_beta <- S_dev_beta_beta(z, cov, cov_j, cov_j2, R, beta, gamma, lambda)
  SR_dev_beta_beta[R == Inf] <- 0
  
  q_dev_b_b <- (1 - c) * (- (1 / (SL - SR)^2) * (SL_dev_beta - SR_dev_beta) * (SL_dev_beta2 - SR_dev_beta2) +
                            1 / (SL - SR) * (SL_dev_beta_beta - SR_dev_beta_beta)
  )
  q_dev_b_b[c == 1] <- 0
  q_dev_b_b
}

k1_dev_lambda_lambda <- function(c, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_lambda <- S_dev_lambda_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda_lambda[L == 0] <- 0
  SR_dev_lambda_lambda <- S_dev_lambda_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda_lambda[R == Inf] <- 0
  
  
  q_dev_lambda_lambda <- (1 - c) * (- (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda)^2 +
                                      1 / (SL - SR) * (SL_dev_lambda_lambda - SR_dev_lambda_lambda))
  q_dev_lambda_lambda[c == 1] <- 0
  q_dev_lambda_lambda
  
}

k1_dev_gamma_gamma <- function(c, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_gamma <- S_dev_gamma_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma_gamma[L == 0] <- 0
  SR_dev_gamma_gamma <- S_dev_gamma_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma_gamma[R == Inf] <- 0
  
  q_dev_gamma_gamma <- (1 - c) * (- (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma)^2 +
                                    1 / (SL - SR) * (SL_dev_gamma_gamma - SR_dev_gamma_gamma))
  q_dev_gamma_gamma[c == 1] <- 0
  q_dev_gamma_gamma
  
}

k1_dev_lambda_beta <- function(c, z, cov, cov_j, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_lambda_beta <- S_dev_lambda_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_lambda_beta[L == 0] <- 0
  SR_dev_lambda_beta <- S_dev_lambda_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_lambda_beta[R == Inf] <- 0
  
  q_dev_lambda_beta <- (1 - c) * (
    - (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_lambda_beta - SR_dev_lambda_beta)
  ) 
  q_dev_lambda_beta[c == 1] <- 0
  q_dev_lambda_beta
  
}

k1_dev_gamma_beta <- function(c, z, cov, cov_j, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_beta <- S_dev_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_beta <- S_dev_gamma_beta(z, cov, cov_j, L, beta, gamma, lambda)
  SL_dev_gamma_beta[L == 0] <- 0
  SR_dev_gamma_beta <- S_dev_gamma_beta(z, cov, cov_j, R, beta, gamma, lambda)
  SR_dev_gamma_beta[R == Inf] <- 0
  
  q_dev_gamma_beta <- (1 - c) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_gamma_beta - SR_dev_gamma_beta)
  ) 
  q_dev_gamma_beta[c == 1] <- 0
  q_dev_gamma_beta
  
}

k1_dev_lambda_gamma <- function(c, z, cov, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z, cov)
  SR <- SurvF(R, gamma, lambda, beta, z, cov)
  
  SL_dev_gamma <- S_dev_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_lambda <- S_dev_lambda(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_gamma <- S_dev_lambda_gamma(z, cov, L, beta, gamma, lambda)
  SL_dev_lambda_gamma[L == 0] <- 0
  SR_dev_lambda_gamma <- S_dev_lambda_gamma(z, cov, R, beta, gamma, lambda)
  SR_dev_lambda_gamma[R == Inf] <- 0
  
  
  q_dev_lambda_gamma <- (1 - c) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_lambda - SR_dev_lambda) +
      1 / (SL - SR) * (SL_dev_lambda_gamma - SR_dev_lambda_gamma)
  ) 
  q_dev_lambda_gamma[c == 1] <- 0
  q_dev_lambda_gamma
  
}

hessianM <- function(a0, aj, beta, gamma, lambda, L, R, c, k, cov_aj, strata, cov_beta){
  # define k1, k0
  k1 = k == 1 # non mixture part
  k0 = k == 0 # mixture part
  
  naj <- ifelse(is.null(cov_aj), 0, ncol(cov_aj))
  nbeta <- ifelse(is.null(cov_beta), 0, ncol(cov_beta))
  n_gamma_lambda <- ncol(strata)
  # calculate pi
  pi <- piF(a0, aj, cov_aj)
  # calculate mixture part
  SR <- SurvF(R, gamma, lambda, beta, strata, cov_beta)
  mx = pi + (1 - pi) * (1 - SR)
  
  hessian_matrix <- matrix(0, ncol = (1 + naj + nbeta + n_gamma_lambda*2), nrow = (1 + naj + nbeta + n_gamma_lambda*2))
  
  dev_a0_a0 <- dev_a0_lambda <- dev_a0_gamma <- dev_a0_beta <- dev_aj_lambda <- dev_aj_gamma <- dev_aj_beta <- 
    dev_lambda_lambda <- dev_lambda_gamma <- dev_gamma_gamma <-
    dev_lambda_beta <- dev_gamma_beta <- dev_a0_aj <-
    dev_beta_beta <- dev_aj_aj <- c()
  
  # for k1, they are all 0. for k0, they are not since mixture part
  dev_a0_gamma[k1] <- dev_a0_lambda[k1] <- dev_a0_beta[k1] <- dev_aj_gamma[k1] <- dev_aj_lambda[k1] <- dev_aj_beta[k1]  <- 0
  
  dev_a0_a0[k1] <- k1_dev_a0_a0(pi, c)[k1]
  dev_a0_a0[k0] <- k0_dev_a0_a0(mx, pi, beta, gamma, lambda, R, strata, cov_beta)[k0]
  
  dev_gamma_gamma[k1] <- k1_dev_gamma_gamma(c, strata, cov_beta, L, R, beta, gamma, lambda)[k1]
  dev_gamma_gamma[k0] <- k0_dev_gamma_gamma(mx, pi, beta, gamma, lambda, R, strata, cov_beta)[k0]
  
  dev_lambda_lambda[k1] <- k1_dev_lambda_lambda(c, strata, cov_beta, L, R, beta, gamma, lambda)[k1]
  dev_lambda_lambda[k0] <- k0_dev_lambda_lambda(mx, pi, beta, gamma, lambda, R, strata, cov_beta)[k0]
  
  dev_lambda_gamma[k1] <- k1_dev_lambda_gamma(c, strata, cov_beta, L, R, beta, gamma, lambda)[k1]
  dev_lambda_gamma[k0] <- k0_dev_lambda_gamma(mx, pi, beta, gamma, lambda, R, strata, cov_beta)[k0]
  
  dev_a0_gamma[k0] <- k0_dev_a0_gamma(mx, pi, beta, gamma, lambda, R, strata, cov_beta)[k0]
  
  dev_a0_lambda[k0] <- k0_dev_a0_lambda(mx, pi, beta, gamma, lambda, R, strata, cov_beta)[k0]
  
  for (j in 1: max(naj, nbeta)) {
    if (j <= naj){
      dev_a0_aj[k1] <- k1_dev_a0_aj(pi, c, cov_aj[, j])[k1]
      dev_a0_aj[k0] <- k0_dev_a0_aj(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_aj[, j])[k0]
      hessian_matrix[1, (1+j)] <- hessian_matrix[(1+j), 1] <- sum(dev_a0_aj)
      
      dev_aj_gamma[k0] <- k0_dev_aj_gamma(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_aj[, j])[k0]
      dev_aj_lambda[k0] <- k0_dev_aj_lambda(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_aj[, j])[k0]
      
      for (q in 1:n_gamma_lambda) {
        hessian_matrix[(1+naj+nbeta+q), (1+j)] <- hessian_matrix[(1+j), (1+naj+nbeta+q)] <- sum(dev_aj_gamma[strata[, q] == 1])
        hessian_matrix[(1+naj+nbeta+n_gamma_lambda+q), (1+j)] <- hessian_matrix[(1+j), (1+naj+nbeta+n_gamma_lambda+q)] <- sum(dev_aj_lambda[strata[, q] == 1])
      }
    }
    
    if (j <= nbeta){
      dev_gamma_beta[k1] <- k1_dev_gamma_beta(c, strata, cov_beta, cov_beta[, j], L, R, beta, gamma, lambda)[k1]
      dev_gamma_beta[k0] <- k0_dev_gamma_beta(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_beta[, j])[k0]
      
      dev_lambda_beta[k1] <- k1_dev_lambda_beta(c, strata, cov_beta, cov_beta[, j], L, R, beta, gamma, lambda)[k1]
      dev_lambda_beta[k0] <- k0_dev_lambda_beta(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_beta[, j])[k0]
      
      for (q in 1:n_gamma_lambda) {
        hessian_matrix[(1+naj+nbeta+q), (1+naj+j)] <- hessian_matrix[(1+naj+j), (1+naj+nbeta+q)] <- sum(dev_gamma_beta[strata[, q] == 1])
        hessian_matrix[(1+naj+nbeta+n_gamma_lambda+q), (1+naj+j)] <- hessian_matrix[(1+naj+j), (1+naj+nbeta+n_gamma_lambda+q)] <- sum(dev_lambda_beta[strata[, q] == 1])
      }
      
      dev_a0_beta[k0] <- k0_dev_a0_beta(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_beta[, j])[k0]
      hessian_matrix[1, (1+naj+j)] <- hessian_matrix[(1+naj+j), 1] <- sum(dev_a0_beta)
    }
    
    for (j2 in 1:max(naj, nbeta)){
      if (j <= naj & j2 <= naj){
        dev_aj_aj[k1] <- k1_dev_aj_aj(pi, c, cov_aj[, j], cov_aj[, j2])[k1]
        dev_aj_aj[k0] <- k0_dev_aj_aj(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_aj[, j], cov_aj[, j2])[k0]
        hessian_matrix[(1+j), (1+j2)] <- hessian_matrix[(1+j2), (1+j)] <- sum(dev_aj_aj)
        
      }
      if (j <= nbeta & j2 <= nbeta){
        dev_beta_beta[k1] <- k1_dev_beta_beta(c, strata, cov_beta, cov_beta[, j], cov_beta[, j2], L, R, beta, gamma, lambda)[k1]
        dev_beta_beta[k0] <- k0_dev_beta_beta(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_beta[, j], cov_beta[, j2])[k0]
        hessian_matrix[(1+naj+j), (1+naj+j2)] <- hessian_matrix[(1+naj+j2), (1+naj+j)] <- sum(dev_beta_beta)
        
      }
      if (j <= naj & j2 <= nbeta){
        dev_aj_beta[k0] <- k0_dev_aj_beta(mx, pi, beta, gamma, lambda, R, strata, cov_beta, cov_beta[, j2], cov_aj[, j])[k0]
        hessian_matrix[(1+j), (1+naj+j2)] <- hessian_matrix[(1+naj+j2), {1+j}] <- sum(dev_aj_beta)
      }
      
    }
  }
  
  
  hessian_matrix[1, 1] <- sum(dev_a0_a0)
  
  for (q in 1:n_gamma_lambda) {
    hessian_matrix[1, (1+naj+nbeta+q)] <- hessian_matrix[(1+naj+nbeta+q), 1] <- sum(dev_a0_gamma[strata[, q] == 1])
    hessian_matrix[1, (1+naj+nbeta+n_gamma_lambda+q)] <- hessian_matrix[(1+naj+nbeta+n_gamma_lambda+q), 1] <- sum(dev_a0_lambda[strata[, q] == 1])
    
    hessian_matrix[(1+naj+nbeta+q), (1+naj+nbeta+q)] <- sum(dev_gamma_gamma[strata[, q] == 1])
    hessian_matrix[(1+naj+nbeta+n_gamma_lambda+q), (1+naj+nbeta+n_gamma_lambda+q)] <- sum(dev_lambda_lambda[strata[, q] == 1])
    hessian_matrix[(1+naj+nbeta+q), (1+naj+nbeta+n_gamma_lambda+q)] <- hessian_matrix[(1+naj+nbeta+n_gamma_lambda+q), (1+naj+nbeta+q)] <- sum(dev_lambda_gamma[strata[, q] == 1])
  }
  
  covMx <- solve(-hessian_matrix) # In = -hessian_matrix, I1 = -hessian_matrix/N
  covMx
  
}
