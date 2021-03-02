library(survival)

beta = 0.7
n=2500

set.seed(20)
z <-  sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
c <-  sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
cov <- sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5)) # binary variable, like gender, cov=0 as female, cov=1 as male
# eventtime <- rweibull(n, shape=lambda*exp(beta*z), scale=gamma)
U <- runif(n) # Generate a random variable U, in  order to find mixture T

##################### Inverse CDF of Weibull ###################
# Inverse CDF of weibull
piF <- function(a0, a1, z){
  linear <- a0 + a1*z
  pi <- exp(linear) / (1 + exp(linear))
  return(pi)
}


eventtimeFc <- function(n, pi, gamma0, gamma1, lambda0, lambda1, beta, cov){
  w <- runif(n)
  eventtime <- c()
    for (i in 1:n) {#print(i)
      if(U[i]<=pi[i]){
        eventtime[i] <- 0 # P(T=0)=pi
      }else{ # P(T>0)=1-pi
        eventtime[i] <- ifelse(z[i] == 0, 
                               (-log(1-w[i])/(lambda0*exp(beta*cov[i])))^(1/gamma0),
                               (-log(1-w[i])/(lambda1*exp(beta*cov[i])))^(1/gamma1)) 
      }
    }
  return(eventtime)
}



################################ pi=0.1 failure_rate=0.9 #######################
pi0 = 0.1
failure_rate = 0.9

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}

a0 = -0.9
a1 = -4
#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta)
        prob <- mean(eventtime<10)
        # if(prob == failure_rate & median(eventtime)<10) {
        if(round(prob, 1) == failure_rate) {
          print(c(g0, g1, l0, l1, max(eventtime)))
        }
      }
    }
  }
  
}

lambda0 = 2
lambda1 = 1.2


################################ pi=0.2 failure_rate=0.9 #######################
pi0 = 0.2
failure_rate = 0.9

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = -0.6
a1 = -3

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

lambda0 = 2.3
lambda1 = 1.1


################################ pi=0.4 failure_rate=0.9 #######################
pi0 = 0.4
failure_rate = 0.9

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = 2.2
a1 = -7

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

lambda0 = 1.2
lambda1 = 0.9

################################ pi=0.1 failure_rate=0.5 #######################
pi0 = 0.1
failure_rate = 0.5

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = -0.9
a1 = -4

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

gamma0 = 0.1
gamma1 = 0.2
lambda0 = 1.0
lambda1 = 0.1

################################ pi=0.2 failure_rate=0.5 #######################
pi0 = 0.2
failure_rate = 0.5

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = -0.1
a1 = -4

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

lambda0 = 0.8
lambda1 = 0.1

################################ pi=0.4 failure_rate=0.5 #######################
pi0 = 0.4
failure_rate = 0.5

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = 1.5
a1 = -4

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.2, 1, by = 0.1)
shapeSeq1 <- seq(0.1, 1, by = 0.1)
scaleSeq0 <- seq(0.2, 1, by = 0.1)
scaleSeq1 <- seq(0.1, 1, by = 0.1)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta, cov)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

lambda0 = 2.000000e-01
lambda1 = 1.000000e-01


################################ pi=0.1 failure_rate=0.3 #######################
pi0 = 0.1
failure_rate = 0.3

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = -0.9
a1 = -4

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta, cov)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

lambda0 = 0.1
lambda1 = 0.1
pi <- piF(a0, a1, z)
eventtime <- eventtimeFc(n, pi, gamma0, gamma1, lambda0, lambda1, beta, cov)
hist(eventtime)
mean(eventtime<10, na.rm = T)

histF(pi, gamma, lambda)

################################ pi=0.2 failure_rate=0.3 #######################
pi0 = 0.2
failure_rate = 0.3

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}


a0 = -0.1
a1 = -4

#  Choose parameters such that 90% of infants  T<10
shapeSeq0 <- seq(0.1, 10, length.out = 100)
shapeSeq1 <- seq(0.1, 10, length.out = 100)
scaleSeq0 <- seq(0.1, 10, length.out = 100)
scaleSeq1 <- seq(0.1, 10, length.out = 100)
pi <- piF(a0, a1, z)

# inverse cdf
for (g0 in shapeSeq0) {
  for (g1 in shapeSeq1) {
    for (l0 in scaleSeq0) {
      for (l1 in scaleSeq1) {
        if (g0 != g1){
          eventtime <- eventtimeFc(n, pi, g0, g1, l0, l1, beta, cov)
          prob <- mean(eventtime<10)
          # if(prob == failure_rate & median(eventtime)<10) {
          if(round(prob, 1) == failure_rate) {
            print(c(g0, g1, l0, l1, max(eventtime)))
          }
        }
      }
    }
  }
  
}

lambda0 = 0.1
lambda1 = 0.1
pi <- piF(a0, a1, z)
eventtime <- eventtimeFc(n, pi, gamma0, gamma1, lambda0, lambda1, beta, cov)
hist(eventtime)
mean(eventtime<10, na.rm = T)

