exp_z_beta <- function(beta, z){
  if(any(is.na(beta))){
    1
  }else{
    exp(z %*% beta)
  }
}

######################################## Hessian matrix ######################################## 
## logistic model ####
pi_dev_a0 <- function(pi){pi - pi^2}

pi_dev_aj <- function(pi, z_aj){z_aj * (pi - pi^2)}

pi_dev_a0_a0 <- function(pi){(1 - 2*pi) * pi_dev_a0(pi)}

pi_dev_aj_aj <- function(pi, z_aj, z_aj2){z_aj * (1 - 2*pi) * pi_dev_aj(pi, z_aj2)}

pi_dev_a0_aj <- function(pi, z_aj){(1 - 2*pi) * pi_dev_aj(pi, z_aj)}

## survival model ####
S_dev_lambda <- function(z, t, beta, gamma, lambda){
  SurvF(t, gamma, lambda, beta, z) * (-t^gamma * exp_z_beta(beta, z))
}

S_dev_gamma <- function(z, t, beta, gamma, lambda){
  SurvF(t, gamma, lambda, beta, z) * (-lambda*exp_z_beta(beta, z)) * t^gamma * log(t)
}

S_dev_beta <- function(z, zj, t, beta, gamma, lambda){
  SurvF(t, gamma, lambda, beta, z) * (-lambda*t^gamma) * exp_z_beta(beta, z) * zj
}

S_dev_lambda_lambda <- function(z, t, beta, gamma, lambda){
  -t^gamma * exp_z_beta(beta, z) * S_dev_lambda(z, t, beta, gamma, lambda)
}

S_dev_gamma_gamma <- function(z, t, beta, gamma, lambda){
  -lambda * exp_z_beta(beta, z) * log(t) * t^gamma *
    (S_dev_gamma(z, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z) * log(t))
}

S_dev_lambda_gamma <- function(z, t, beta, gamma, lambda){
  -exp_z_beta(beta, z) * t^gamma * (S_dev_gamma(z, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z)*log(t))
}

S_dev_lambda_beta <- function(z, zj, t, beta, gamma, lambda){
  -t^gamma * exp_z_beta(beta, z) * (S_dev_beta(z, zj, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z)*zj)
}

S_dev_gamma_beta <- function(z, zj, t, beta, gamma, lambda){
  -lambda * t^gamma * log(t) * exp_z_beta(beta, z) * (S_dev_beta(z, zj, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z)*zj)
}

S_dev_beta_beta <- function(z, zj, zj2, t, beta, gamma, lambda){
  -lambda * t^gamma * zj * exp_z_beta(beta, z) * (S_dev_beta(z, zj2, t, beta, gamma, lambda) + SurvF(t, gamma, lambda, beta, z)*zj2)
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

Q_dev_beta <- function(Ei, z, zj, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  q_dev_b <- (1 - Ei) * (1 / (SL - SR)) * (SL_dev_beta - SR_dev_beta)
  q_dev_b[Ei == 1] <- 0
  q_dev_b
}

Q_dev_lambda <- function(Ei, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  q_dev_lambda <- (1 - Ei) * (1 / (SL - SR)) * (SL_dev_lambda - SR_dev_lambda)
  q_dev_lambda[Ei == 1] <- 0
  q_dev_lambda
  
}

Q_dev_gamma <- function(Ei, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  q_dev_gamma <- (1 - Ei) * (1 / (SL - SR)) * (SL_dev_gamma - SR_dev_gamma)
  q_dev_gamma[Ei == 1] <- 0
  q_dev_gamma
}

Q_dev_beta_beta <- function(Ei, z, zj, zj2, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_beta2 <- S_dev_beta(z, zj2, L, beta, gamma, lambda)
  SL_dev_beta2[L == 0] <- 0
  SR_dev_beta2 <- S_dev_beta(z, zj2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  
  SL_dev_beta_beta <- S_dev_beta_beta(z, zj, zj2, L, beta, gamma, lambda)
  SL_dev_beta_beta[L == 0] <- 0
  SR_dev_beta_beta <- S_dev_beta_beta(z, zj, zj2, R, beta, gamma, lambda)
  SR_dev_beta_beta[R == Inf] <- 0
  
  q_dev_b_b <- (1 - Ei) * (- (1 / (SL - SR)^2) * (SL_dev_beta - SR_dev_beta) * (SL_dev_beta2 - SR_dev_beta2) +
                             1 / (SL - SR) * (SL_dev_beta_beta - SR_dev_beta_beta)
  )
  q_dev_b_b[Ei == 1] <- 0
  q_dev_b_b
}

Q_dev_lambda_lambda <- function(Ei, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_lambda <- S_dev_lambda_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda_lambda[L == 0] <- 0
  SR_dev_lambda_lambda <- S_dev_lambda_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda_lambda[R == Inf] <- 0
  
  
  q_dev_lambda_lambda <- (1 - Ei) * (- (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda)^2 +
                                       1 / (SL - SR) * (SL_dev_lambda_lambda - SR_dev_lambda_lambda))
  q_dev_lambda_lambda[Ei == 1] <- 0
  q_dev_lambda_lambda
  
}

Q_dev_gamma_gamma <- function(Ei, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_gamma <- S_dev_gamma_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma_gamma[L == 0] <- 0
  SR_dev_gamma_gamma <- S_dev_gamma_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma_gamma[R == Inf] <- 0
  
  q_dev_gamma_gamma <- (1 - Ei) * (- (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma)^2 +
                                     1 / (SL - SR) * (SL_dev_gamma_gamma - SR_dev_gamma_gamma))
  q_dev_gamma_gamma[Ei == 1] <- 0
  q_dev_gamma_gamma
  
}

Q_dev_lambda_beta <- function(Ei, z, zj, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_lambda_beta <- S_dev_lambda_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_lambda_beta[L == 0] <- 0
  SR_dev_lambda_beta <- S_dev_lambda_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_lambda_beta[R == Inf] <- 0
  
  q_dev_lambda_beta <- (1 - Ei) * (
    - (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_lambda_beta - SR_dev_lambda_beta)
  ) 
  q_dev_lambda_beta[Ei == 1] <- 0
  q_dev_lambda_beta
  
}

Q_dev_gamma_beta <- function(Ei, z, zj, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_beta <- S_dev_gamma_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_gamma_beta[L == 0] <- 0
  SR_dev_gamma_beta <- S_dev_gamma_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_gamma_beta[R == Inf] <- 0
  
  q_dev_gamma_beta <- (1 - Ei) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_gamma_beta - SR_dev_gamma_beta)
  ) 
  q_dev_gamma_beta[Ei == 1] <- 0
  q_dev_gamma_beta
  
}

Q_dev_lambda_gamma <- function(Ei, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_gamma <- S_dev_lambda_gamma(z, L, beta, gamma, lambda)
  SL_dev_lambda_gamma[L == 0] <- 0
  SR_dev_lambda_gamma <- S_dev_lambda_gamma(z, R, beta, gamma, lambda)
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
k0_dev_a0_a0 <- function(mx, pi, beta, gamma, lambda, R, z){
  SR <- SurvF(R, gamma, lambda, beta, z)
  (1 / mx) * SR * (pi_dev_a0_a0(pi) - (1 / mx) * SR * pi_dev_a0(pi)^2)
}

k0_dev_aj_aj <- function(mx, pi, beta, gamma, lambda, R, z, z_aj, z_aj2){
  SR <- SurvF(R, gamma, lambda, beta, z)
  SR * (1 / mx) * (pi_dev_aj_aj(pi, z_aj, z_aj2) - (1 / mx) * SR * pi_dev_aj(pi, z_aj) * pi_dev_aj(pi, z_aj2))
}

k0_dev_a0_aj <- function(mx, pi, beta, gamma, lambda, R, z, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z)
  SR * (1/mx) * (pi_dev_a0_aj(pi, z_aj) - (1/mx) * SR * pi_dev_a0(pi) * pi_dev_aj(pi, z_aj))
}

k0_dev_beta_beta <- function(mx, pi, beta, gamma, lambda, R, z, zj, zj2){
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  SR_dev_beta2 <- S_dev_beta(z, zj2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  SR_dev_beta_beta <- S_dev_beta_beta(z, zj, zj2, R, beta, gamma, lambda)
  SR_dev_beta_beta[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx)*(1-pi) * SR_dev_beta * SR_dev_beta2 + SR_dev_beta_beta)
}

k0_dev_gamma_gamma <- function(mx, pi, beta, gamma, lambda, R, z){
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  SR_dev_gamma_gamma <- S_dev_gamma_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma_gamma[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_gamma^2 + SR_dev_gamma_gamma)
}

k0_dev_lambda_lambda <- function(mx, pi, beta, gamma, lambda, R, z){
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  SR_dev_lambda_lambda <- S_dev_lambda_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda_lambda[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_lambda^2 + SR_dev_lambda_lambda)
}

k0_dev_gamma_beta <- function(mx, pi, beta, gamma, lambda, R, z, zj){
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  SR_dev_gamma_beta <- S_dev_gamma_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_gamma_beta[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_gamma * SR_dev_beta + SR_dev_gamma_beta)
}

k0_dev_lambda_beta <- function(mx, pi, beta, gamma, lambda, R, z, zj){
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  SR_dev_lambda_beta <- S_dev_lambda_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_lambda_beta[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_lambda * SR_dev_beta + SR_dev_lambda_beta)
}

k0_dev_lambda_gamma <- function(mx, pi, beta, gamma, lambda, R, z){
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  SR_dev_lambda_gamma <- S_dev_lambda_gamma(z, R, beta, gamma, lambda)
  SR_dev_lambda_gamma[R == Inf] <- 0
  
  -(1/mx) * (1-pi) * ((1/mx) * (1-pi) * SR_dev_lambda * SR_dev_gamma + SR_dev_lambda_gamma)
}

k0_dev_a0_beta <- function(mx, pi, beta, gamma, lambda, R, z, zj){
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  pi_dev_a0(pi) * (1 / mx) * SR_dev_beta * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_a0_gamma <- function(mx, pi, beta, gamma, lambda, R, z){
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  pi_dev_a0(pi) * (1 / mx) * SR_dev_gamma * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_a0_lambda <- function(mx, pi, beta, gamma, lambda, R, z){
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  pi_dev_a0(pi) * (1 / mx) * SR_dev_lambda * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_aj_beta <- function(mx, pi, beta, gamma, lambda, R, z, zj2, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SR_dev_beta2 <- S_dev_beta(z, zj2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  
  pi_dev_aj(pi, z_aj) * (1 / mx) * SR_dev_beta2 * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_aj_gamma <- function(mx, pi, beta, gamma, lambda, R, z, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  pi_dev_aj(pi, z_aj) * (1 / mx) * SR_dev_gamma * (1 + (1 / mx) * (1 - pi) * SR)
}

k0_dev_aj_lambda <- function(mx, pi, beta, gamma, lambda, R, z, z_aj){
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
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

k1_dev_beta_beta <- function(c, z, zj, zj2, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_beta2 <- S_dev_beta(z, zj2, L, beta, gamma, lambda)
  SL_dev_beta2[L == 0] <- 0
  SR_dev_beta2 <- S_dev_beta(z, zj2, R, beta, gamma, lambda)
  SR_dev_beta2[R == Inf] <- 0
  
  SL_dev_beta_beta <- S_dev_beta_beta(z, zj, zj2, L, beta, gamma, lambda)
  SL_dev_beta_beta[L == 0] <- 0
  SR_dev_beta_beta <- S_dev_beta_beta(z, zj, zj2, R, beta, gamma, lambda)
  SR_dev_beta_beta[R == Inf] <- 0
  
  q_dev_b_b <- (1 - c) * (- (1 / (SL - SR)^2) * (SL_dev_beta - SR_dev_beta) * (SL_dev_beta2 - SR_dev_beta2) +
                            1 / (SL - SR) * (SL_dev_beta_beta - SR_dev_beta_beta)
  )
  q_dev_b_b[c == 1] <- 0
  q_dev_b_b
}

k1_dev_lambda_lambda <- function(c, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_lambda <- S_dev_lambda_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda_lambda[L == 0] <- 0
  SR_dev_lambda_lambda <- S_dev_lambda_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda_lambda[R == Inf] <- 0
  
  
  q_dev_lambda_lambda <- (1 - c) * (- (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda)^2 +
                                      1 / (SL - SR) * (SL_dev_lambda_lambda - SR_dev_lambda_lambda))
  q_dev_lambda_lambda[c == 1] <- 0
  q_dev_lambda_lambda
  
}

k1_dev_gamma_gamma <- function(c, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_gamma <- S_dev_gamma_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma_gamma[L == 0] <- 0
  SR_dev_gamma_gamma <- S_dev_gamma_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma_gamma[R == Inf] <- 0
  
  q_dev_gamma_gamma <- (1 - c) * (- (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma)^2 +
                                    1 / (SL - SR) * (SL_dev_gamma_gamma - SR_dev_gamma_gamma))
  q_dev_gamma_gamma[c == 1] <- 0
  q_dev_gamma_gamma
  
}

k1_dev_lambda_beta <- function(c, z, zj, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_lambda_beta <- S_dev_lambda_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_lambda_beta[L == 0] <- 0
  SR_dev_lambda_beta <- S_dev_lambda_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_lambda_beta[R == Inf] <- 0
  
  q_dev_lambda_beta <- (1 - c) * (
    - (1 / (SL - SR)^2) * (SL_dev_lambda - SR_dev_lambda) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_lambda_beta - SR_dev_lambda_beta)
  ) 
  q_dev_lambda_beta[c == 1] <- 0
  q_dev_lambda_beta
  
}

k1_dev_gamma_beta <- function(c, z, zj, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_beta <- S_dev_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_beta[L == 0] <- 0
  SR_dev_beta <- S_dev_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_beta[R == Inf] <- 0
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_gamma_beta <- S_dev_gamma_beta(z, zj, L, beta, gamma, lambda)
  SL_dev_gamma_beta[L == 0] <- 0
  SR_dev_gamma_beta <- S_dev_gamma_beta(z, zj, R, beta, gamma, lambda)
  SR_dev_gamma_beta[R == Inf] <- 0
  
  q_dev_gamma_beta <- (1 - c) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_beta - SR_dev_beta) +
      1 / (SL - SR) * (SL_dev_gamma_beta - SR_dev_gamma_beta)
  ) 
  q_dev_gamma_beta[c == 1] <- 0
  q_dev_gamma_beta
  
}

k1_dev_lambda_gamma <- function(c, z, L, R, beta, gamma, lambda){
  SL <- SurvF(L, gamma, lambda, beta, z)
  SR <- SurvF(R, gamma, lambda, beta, z)
  
  SL_dev_gamma <- S_dev_gamma(z, L, beta, gamma, lambda)
  SL_dev_gamma[L == 0] <- 0
  SR_dev_gamma <- S_dev_gamma(z, R, beta, gamma, lambda)
  SR_dev_gamma[R == Inf] <- 0
  
  SL_dev_lambda <- S_dev_lambda(z, L, beta, gamma, lambda)
  SL_dev_lambda[L == 0] <- 0
  SR_dev_lambda <- S_dev_lambda(z, R, beta, gamma, lambda)
  SR_dev_lambda[R == Inf] <- 0
  
  SL_dev_lambda_gamma <- S_dev_lambda_gamma(z, L, beta, gamma, lambda)
  SL_dev_lambda_gamma[L == 0] <- 0
  SR_dev_lambda_gamma <- S_dev_lambda_gamma(z, R, beta, gamma, lambda)
  SR_dev_lambda_gamma[R == Inf] <- 0
  
  
  q_dev_lambda_gamma <- (1 - c) * (
    - (1 / (SL - SR)^2) * (SL_dev_gamma - SR_dev_gamma) * 
      (SL_dev_lambda - SR_dev_lambda) +
      1 / (SL - SR) * (SL_dev_lambda_gamma - SR_dev_lambda_gamma)
  ) 
  q_dev_lambda_gamma[c == 1] <- 0
  q_dev_lambda_gamma
  
}

hessianEM <- function(a0, aj, beta, gamma, lambda, L, R, c, k, z_aj, z){
  # define k1, k0
  k1 = k == 1 # non mixture part
  k0 = k == 0 # mixture part
  
  naj <- ifelse(is.null(z_aj), 0, ncol(z_aj))
  nbeta <- ifelse(is.null(z), 0, ncol(z))
  # calculate mixture component Ei
  # Ei <- EiF(beta, gamma, lambda, a0, aj, L, R, k, c, z_aj, z)
  # calculate pi
  pi <- piF(a0, aj, z_aj)
  # calculate mixture part
  SR <- SurvF(R, gamma, lambda, beta, z)
  mx = pi + (1 - pi) * (1 - SR)
  
  hessian_matrix <- matrix(0, ncol = (1 + naj + nbeta + 1 + 1), nrow = (1 + naj + nbeta + 1 + 1))
  
  dev_a0_a0 <- dev_a0_lambda <- dev_a0_gamma <- dev_a0_beta <- dev_aj_lambda <- dev_aj_gamma <- dev_aj_beta <- 
    dev_lambda_lambda <- dev_lambda_gamma <- dev_gamma_gamma <-
    dev_lambda_beta <- dev_gamma_beta <- dev_a0_aj <-
    dev_beta_beta <- dev_aj_aj <- c()
  
  # for k1, they are all 0. for k0, they are not since mixture part
  dev_a0_gamma[k1] <- dev_a0_lambda[k1] <- dev_a0_beta[k1] <- dev_aj_gamma[k1] <- dev_aj_lambda[k1] <- dev_aj_beta[k1]  <- 0
  
  dev_a0_a0[k1] <- k1_dev_a0_a0(pi, c)[k1]
  dev_a0_a0[k0] <- k0_dev_a0_a0(mx, pi, beta, gamma, lambda, R, z)[k0]
  
  dev_gamma_gamma[k1] <- k1_dev_gamma_gamma(c, z, L, R, beta, gamma, lambda)[k1]
  dev_gamma_gamma[k0] <- k0_dev_gamma_gamma(mx, pi, beta, gamma, lambda, R, z)[k0]
  
  dev_lambda_lambda[k1] <- k1_dev_lambda_lambda(c, z, L, R, beta, gamma, lambda)[k1]
  dev_lambda_lambda[k0] <- k0_dev_lambda_lambda(mx, pi, beta, gamma, lambda, R, z)[k0]
  
  dev_lambda_gamma[k1] <- k1_dev_lambda_gamma(c, z, L, R, beta, gamma, lambda)[k1]
  dev_lambda_gamma[k0] <- k0_dev_lambda_gamma(mx, pi, beta, gamma, lambda, R, z)[k0]
  
  dev_a0_gamma[k0] <- k0_dev_a0_gamma(mx, pi, beta, gamma, lambda, R, z)[k0]
  
  dev_a0_lambda[k0] <- k0_dev_a0_lambda(mx, pi, beta, gamma, lambda, R, z)[k0]
  
  if (naj > 0 | nbeta > 0){
    for (j in 1: max(naj, nbeta)) {
      if (naj > 0 & j <= naj){
        dev_a0_aj[k1] <- k1_dev_a0_aj(pi, c, z_aj[, j])[k1]
        dev_a0_aj[k0] <- k0_dev_a0_aj(mx, pi, beta, gamma, lambda, R, z, z_aj[, j])[k0]
        hessian_matrix[1, (1+j)] <- hessian_matrix[(1+j), 1] <- sum(dev_a0_aj)
        
        dev_aj_gamma[k0] <- k0_dev_aj_gamma(mx, pi, beta, gamma, lambda, R, z, z_aj[, j])[k0]
        hessian_matrix[(1+naj+nbeta+1), (1+j)] <- hessian_matrix[(1+j), (1+naj+nbeta+1)] <- sum(dev_aj_gamma)
        
        dev_aj_lambda[k0] <- k0_dev_aj_lambda(mx, pi, beta, gamma, lambda, R, z, z_aj[, j])[k0]
        hessian_matrix[(1+naj+nbeta+2), (1+j)] <- hessian_matrix[(1+j), (1+naj+nbeta+2)] <- sum(dev_aj_lambda)
        
      }
      
      if (nbeta > 0 & j <= nbeta){
        dev_gamma_beta[k1] <- k1_dev_gamma_beta(c, z, z[, j], L, R, beta, gamma, lambda)[k1]
        dev_gamma_beta[k0] <- k0_dev_gamma_beta(mx, pi, beta, gamma, lambda, R, z, z[, j])[k0]
        hessian_matrix[(1+naj+nbeta+1), (1+naj+j)] <- hessian_matrix[(1+naj+j), (1+naj+nbeta+1)] <- sum(dev_gamma_beta)
        
        dev_lambda_beta[k1] <- k1_dev_lambda_beta(c, z, z[, j], L, R, beta, gamma, lambda)[k1]
        dev_lambda_beta[k0] <- k0_dev_lambda_beta(mx, pi, beta, gamma, lambda, R, z, z[, j])[k0]
        hessian_matrix[(1+naj+nbeta+2), (1+naj+j)] <- hessian_matrix[(1+naj+j), (1+naj+nbeta+2)] <- sum(dev_lambda_beta)
        
        dev_a0_beta[k0] <- k0_dev_a0_beta(mx, pi, beta, gamma, lambda, R, z, z[, j])[k0]
        hessian_matrix[1, (1+naj+j)] <- hessian_matrix[(1+naj+j), 1] <- sum(dev_a0_beta)
      }
      
      for (j2 in 1:max(naj, nbeta)){
        if (naj > 0 & j <= naj & j2 <= naj){
          dev_aj_aj[k1] <- k1_dev_aj_aj(pi, c, z_aj[, j], z_aj[, j2])[k1]
          dev_aj_aj[k0] <- k0_dev_aj_aj(mx, pi, beta, gamma, lambda, R, z, z_aj[, j], z_aj[, j2])[k0]
          hessian_matrix[(1+j), (1+j2)] <- hessian_matrix[(1+j2), (1+j)] <- sum(dev_aj_aj)
          
        }
        if (nbeta > 0 & j <= nbeta & j2 <= nbeta){
          dev_beta_beta[k1] <- k1_dev_beta_beta(c, z, z[, j], z[, j2], L, R, beta, gamma, lambda)[k1]
          dev_beta_beta[k0] <- k0_dev_beta_beta(mx, pi, beta, gamma, lambda, R, z, z[, j], z[, j2])[k0]
          hessian_matrix[(1+naj+j), (1+naj+j2)] <- hessian_matrix[(1+naj+j2), (1+naj+j)] <- sum(dev_beta_beta)
          
        }
        if (naj > 0 & nbeta > 0 & j <= naj & j2 <= nbeta){
          dev_aj_beta[k0] <- k0_dev_aj_beta(mx, pi, beta, gamma, lambda, R, z, z[, j2], z_aj[, j])[k0]
          hessian_matrix[(1+j), (1+naj+j2)] <- hessian_matrix[(1+naj+j2), {1+j}] <- sum(dev_aj_beta)
        }
        
      }
    }
    
  }
  
  hessian_matrix[1, 1] <- sum(dev_a0_a0)
  hessian_matrix[(1+naj+nbeta+1), (1+naj+nbeta+1)] <- sum(dev_gamma_gamma)
  hessian_matrix[(1+naj+nbeta+2), (1+naj+nbeta+2)] <- sum(dev_lambda_lambda)
  hessian_matrix[(1+naj+nbeta+1), (1+naj+nbeta+2)] <- hessian_matrix[(1+naj+nbeta+2), (1+naj+nbeta+1)] <- sum(dev_lambda_gamma)
  
  covMx <- solve(-hessian_matrix)
  covMx
  
}


