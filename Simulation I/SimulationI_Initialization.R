library(icensmis)
library(straweib)
library(survival)


n=2500

set.seed(20)
z <-  sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
c <-  sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
U <- runif(n) # Generate a random variable U, in  order to find mixture T

##################### Inverse CDF of Weibull ###################
# Inverse CDF of weibull
piF <- function(a0, a1, z){
  linear <- a0 + a1*z
  pi <- exp(linear) / (1 + exp(linear))
  return(pi)
}

# 
eventtimeFc <- function(n, pi, gamma, lambda){
  w <- runif(n)
  eventtime <- c()
    for (i in 1:n) {#print(i)
      if(U[i]<=pi[i]){
        eventtime[i] <- 0 # P(T=0)=pi
      }else{ # P(T>0)=1-pi
        eventtime[i] <- (-log(1-w[i])/(lambda*exp(beta*z[i])))^(1/gamma)
      }
    }
  return(eventtime)
}

histF <- function(pi, gamma, lambda){
  beta_hat <- replicate(1000, {
    eventtime <- eventtimeFc(n, pi, gamma, lambda)
    tau <- c(1,4,10)
    data <- data.frame(ID = rep(1:n, each=length(tau)), 
                       testtime = rep(tau, n),  
                       eventtime = rep(eventtime, each=length(tau)),
                       z = rep(z, each=length(tau)),
                       c = rep(c, each=length(tau)))
    data$result <- ifelse(data$eventtime >  data$testtime, 0, 1)
    id <- unique(data$ID)
    L0 <-
      sapply(id, function(x)
        tail(data[data$ID == x &
                    data$result == 0, ], 1)$testtime) # choose the last 0
    L <-
      sapply(1:length(id), function(x)
        ifelse(length(L0[[x]]) == 0, 0, L0[[x]]))
    L[data$eventtime != 0 & data$c == 0] <- 0
    R0 <-
      sapply(id, function(x)
        head(data[data$ID == x &
                    data$result == 1, ], 1)$testtime) # choose the first 1
    R <-
      sapply(1:length(id), function(x)
        ifelse(length(R0[[x]]) == 0, Inf, R0[[x]]))
    
    data2 <- data.frame(id = id, 
                        Left = L, 
                        Right = R,
                        z = z)
    data2$Left <- ifelse(is.na(data2$Left), 0, data2$Left)
    data2$Right <- ifelse(is.na(data2$Right), Inf, data2$Right)
    fit <- icweib(Left, Right, data2, strata = 1, ~z)
    as.numeric(fit$coef[1])
  })
  
  label <- paste0("./pic/II_pi_", pi, "_beta.png")
  png(filename=label)
  hist(beta_hat)
  dev.off()
  loc <- paste0("./pic/II_pi_", pi, "_beta.csv")
  write.csv(as.matrix(summary(beta_hat)), file = loc)
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
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.290028
lambda=1.090108


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

a0 = -0.1
a1 = -4
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.550054 
lambda=0.450044

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

a0 = 1.6
a1 = -4.4
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma = 3.000200e-02
lambda = 7.100700e-01

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
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.690068
lambda=0.070006

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
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.570056
lambda=0.070006

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

a0 = 1.7
a1 = -4.5
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    # eventtime <- eventtimeFc(n, pi, k,j)
    eventtime <- eventtimeFc(n, z, pi, k,j, beta=0.7)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma = 2.500240e-01 
lambda = 5.000400e-02

################################ pi=0.05 failure_rate=0.3 #######################
pi0 = 0.05
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
a1 = -4.2
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}
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
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma = 0.530052
lambda=0.050004

################################ pi=0.2 failure_rate=0.3 #######################
pi0 = 0.2
failure_rate = 0.3

a0_seq <- seq(-20, 20, 0.1)
a1_seq <- seq(-20, 20, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}

a0 = -4.2
a1 = 4
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=3.300320e-01
lambda=5.000400e-02

################################ pi=0.4 failure_rate=0.3 (Not Found!) #######################
pi0 = 0.4
failure_rate = 0.3

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 1) == pi0){print(cbind(p, q))}
  }
}

a0 = 2.1
a1 = -8.1
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.450044
lambda=0.190018



################################ pi=0.02 failure_rate=0.1 #######################
pi0 = 0.02
failure_rate = 0.1

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 2) == pi0){print(cbind(p, q))}
  }
}

a0 = -3
a1 = -3.7
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.21002
lambda=0.030002

################################ pi=0.04 failure_rate=0.1 #######################
pi0 = 0.04
failure_rate = 0.1

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 2) == pi0){print(cbind(p, q))}
  }
}

a0 = -2.4
a1 = -3
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.890088 
lambda=0.010000


pi <- piF(a0, a1, z)
eventtime <- eventtimeFc(n, pi, gamma, lambda)
hist(eventtime)
mean(eventtime<10, na.rm = T)

################################ pi=0.05 failure_rate=0.2 #######################
pi0 = 0.05
failure_rate = 0.2

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 2) == pi0){print(cbind(p, q))}
  }
}

a0 = -2.1
a1 = -7
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.710070 
lambda=0.030002


pi <- piF(a0, a1, z)
eventtime <- eventtimeFc(n, pi, gamma, lambda)
hist(eventtime)
mean(eventtime<10, na.rm = T)

################################ pi=0.10 failure_rate=0.2 #######################
pi0 = 0.10
failure_rate = 0.2

a0_seq <- seq(-10, 10, 0.1)
a1_seq <- seq(-10, 10, 0.1)
for (p in a0_seq) {
  for (q in a1_seq){
    pi <- piF(p, q, z)
    if (round(mean(pi), 2) == pi0){print(cbind(p, q))}
  }
}

a0 = -1.4
a1 = -4
#  Choose parameters such that 90% of infants  T<10
shapeSeq <- seq(0.01, 100, length.out = 5000)
scaleSeq <- seq(0.01, 100, length.out = 5000)
pi <- piF(a0, a1, z)

# inverse cdf
for (k in shapeSeq) {
  for (j in scaleSeq) {
    eventtime <- eventtimeFc(n, pi, k,j)
    prob <- mean(eventtime<10)
    # if(prob == failure_rate & median(eventtime)<10) {
    if(round(prob, 1) == failure_rate) {
      print(c(k,j, max(eventtime)))
    }
  }
}

gamma=0.610060 
lambda=0.030002


pi <- piF(a0, a1, z)
eventtime <- eventtimeFc(n, pi, gamma, lambda)
hist(eventtime)
mean(eventtime<10, na.rm = T)

