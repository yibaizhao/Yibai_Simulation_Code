library(ggplot2)
# install.packages("devtools")
# devtools::install_github("teunbrand/ggh4x", force = TRUE)
library(ggh4x)
library(plyr)
require(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
######################################## Scenario I ######################################## 
eventtimeFc <- function(n, pi, gamma, lambda, z){
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

piF <- function(a0, a1, z){
  linear <- a0 + a1*z
  pi <- exp(linear) / (1 + exp(linear))
  return(pi)
}

set.seed(20)
n <- 2500
z <-  sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
U <- runif(n) # Generate a random variable U, in  order to find mixture T

pi0 <- head(rep(c(0.1, 0.2, 0.4), 3), -1)
failure_rate <- fr <- c(rep(0.9, 3), rep(0.5, 3), rep(0.3, 2))
value = matrix(round(c(-0.9, -4, 0.7, 0.290028, 0.910090,
          -0.1, -4, 0.7, 0.550054, 0.450044,
          1.6, -4.4, 0.7, 3.000200e-02, 7.100700e-01,
          -0.9, -4, 0.7, 0.690068, 0.070006,
          -0.1, -4, 0.7, 0.570056, 0.070006,
          1.7, -4.5, 0.7, 2.500240e-01 , 5.000400e-02,
          -0.9, -4, 0.7, 0.530052, 0.050004,
          -4.2, 4, 0.7, 3.300320e-01, 5.000400e-02), 2), byrow = TRUE, ncol = 5)
a0 <- value[,1]
aj <- value[,2]
gamma <- value[,4]
lambda <- value[,5]
beta <- 0.7

simPlots <- function(n, pi0, failure_rate, a0, aj, beta, gamma, lambda, z){
  data_All <- c()
  for (i in 1:length(a0)) {
    pi <- piF(a0[i], aj[i], z)
    t <- eventtimeFc(n, pi, gamma[i], lambda[i], z)
    # setting_title <- paste0('failure rate: ', failure_rate[i],
    #                         ', pi: ', pi0[i],
    #                         '\na0=', a0[i], ', aj=', aj[i], ', beta=0.7', ', gamma=', gamma[i], ', lambda=', lambda[i])
    setting_title <- paste0('a0=', a0[i], ', aj=', aj[i], ', \nbeta=0.7', ', gamma=', gamma[i], ', lambda=', lambda[i])
    data <- data.frame(x = 1:n, t = t)
    data$t_ceiling <- ceiling(data$t)
    data$settings <- setting_title
    data$group <- factor(ifelse(t == 0, 'Test positive at time zero', 'Test positive after time zero'), 
                         levels = c('Test positive at time zero', 'Test positive after time zero'))
    data$z <- ifelse(z == 0, 'Group 0', 'Group 1')
    data$fr <- paste0('Failure rate: ', failure_rate[i])
    data$pi <- factor(paste0('pi: ', pi0[i]), levels = paste0('pi: ', unique(pi0)))
    # data$pi <- factor(paste0(expression(pi), ': ', pi0[i]), levels = paste0(expression(pi), ': ', unique(pi0)))
    data_All <- rbind(data_All, data)
  }
  data_All$settings <- factor(data_All$settings, levels = unique(data_All$settings))
  gList <- list()
  gList[[1]] <- 
    ggplot(subset(data_All[data_All$fr == paste0('Failure rate: ', unique(fr)[1]),], t_ceiling <= 10), aes(x = t_ceiling, fill = group)) +  
      geom_histogram(aes(y = ..count..  / n), bins = 30) +
      scale_x_continuous(breaks = seq(0, 10, 1)) +
      geom_vline(xintercept=10, col = 'red', linetype="dashed") +
      # facet_nested(fr ~ pi + settings + z, scales = 'free', nest_line = TRUE, bleed = TRUE) +
      facet_nested(fr ~ pi + z, scales = 'free', nest_line = TRUE, bleed = TRUE) +
      labs(title = "Distribution of Age (days) at which disease is reported",
           x = '',
           y = '') +
      theme(legend.position="top", legend.title = element_blank())

  gList[[2]] <- 
    ggplot(subset(data_All[data_All$fr == paste0('Failure rate: ', unique(fr)[2]),], t_ceiling <= 10), aes(x = t_ceiling, fill = group)) +  
    geom_histogram(aes(y = ..count..  / n), bins = 30) +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    geom_vline(xintercept=10, col = 'red', linetype="dashed") +
    facet_nested(fr ~ pi + z, scales = 'free', nest_line = TRUE, bleed = TRUE) +
    labs(title = "",
         x = '',
         y = 'Proportion of Positive Tests') +
    theme(legend.position="none")
  
  gList[[3]] <- 
    ggplot(subset(data_All[data_All$fr == paste0('Failure rate: ', unique(fr)[3]),], t_ceiling <= 10), aes(x = t_ceiling, fill = group)) +  
    geom_histogram(aes(y = ..count..  / n), bins = 30) +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    geom_vline(xintercept=10, col = 'red', linetype="dashed") +
    facet_nested(fr ~ pi + z, scales = 'free', nest_line = TRUE, bleed = TRUE) +
    labs(title = "",
         x = 'Age (days) at which test is administered',
         y = '') +
    theme(legend.position="none")  
  gAll <- marrangeGrob(gList, nrow=3, ncol=1, top = NULL)
  gAll
  ggsave("scenarioI_V3.png", plot = gAll, width = 12, height = 8)
}

simPlots(n, pi0, fr, a0, aj, beta, gamma, lambda, z)

######################################## Scenario II ######################################## 
piFII <- function(a0, aj, z){ # z is an indicator of treatment group
  if(is.null(z)){
    linear <- a0 + matrix(0, nrow = nrow(data))
  }else{
    linear <- a0 + z %*% aj
  }
  pi <- exp(linear) / (1 + exp(linear))
  return(pi)
}


eventtimeFcII <- function(n, z, cov, pi, gamma, lambda, beta){ # gamma and lambda are parameters from Weibull(gamma, lambda)
  U <- runif(n) # Generate a uniformly distributed random variable U, in  order to find mixture T
  w <- runif(n) # Generate a uniformly distributed random variable w, representing F, where F=1-S
  eventtime <- c() 
  for (i in 1:n) {
    if(U[i]<=pi[i]){ # if U<=pi, 
      eventtime[i] <- 0 # P(T=0)=pi
    }else{ # if U>pi, P(T>0)=1-pi
      eventtime[i] <- ifelse(z[i] == 0, 
                             (-log(1-w[i])/(lambda[1]*exp(cov[i,] %*% beta)))^(1/gamma[1]),
                             (-log(1-w[i])/(lambda[2]*exp(cov[i,] %*% beta)))^(1/gamma[2])) # Derive eventtime from inverse CDF of Weibull
    }
  }
  return(eventtime)
}

simPlotsII <- function(n, pi0, failure_rate, a0, aj, beta, gamma, lambda, z, cov){
  data_All <- c()
  for (i in 1:length(a0)) {
    pi <- piFII(a0[i], aj[i], z)
    t <- eventtimeFcII(n, z, cov, pi, gamma[i,], lambda[i,], beta)
    setting_title <- paste0('a0=', a0[i], ', aj=', aj[i], ', \nbeta=0.7', 
                            paste0(paste0(', gamma', c(1:length(gamma[i,])), '=', gamma[i,]), collapse = ''),
                            paste0(paste0(', lambda', c(1:length(lambda[i,])), '=', lambda[i,]), collapse = ''))
    data <- data.frame(x = 1:n, t = t)
    data$t_ceiling <- ceiling(data$t)
    data$settings <- setting_title
    data$group <- factor(ifelse(t == 0, 'Test positive at time zero', 'Test positive after time zero'), 
                         levels = c('Test positive at time zero', 'Test positive after time zero'))
    data$z <- ifelse(z == 0, 'Group 0', 'Group 1')
    data$cov <- paste0('cov=', cov)
    data$fr <- paste0('Failure rate: ', failure_rate[i])
    data$pi <- factor(paste0('pi: ', pi0[i]), levels = paste0('pi: ', unique(pi0)))
    data_All <- rbind(data_All, data)
  }
  data_All$settings <- factor(data_All$settings, levels = unique(data_All$settings))

  gList <- list()
  gList[[1]] <- 
    ggplot(subset(data_All[data_All$fr == paste0('Failure rate: ', unique(fr)[1]),], t_ceiling <= 10), aes(x = t_ceiling, fill = group)) +  
    geom_histogram(aes(y = ..count..  / n), bins = 30) +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    geom_vline(xintercept=10, col = 'red', linetype="dashed") +
    facet_nested(fr ~ pi + z + cov, scales = 'free', nest_line = TRUE, bleed = TRUE) +
    labs(title = "Distribution of Age (days) at which disease is reported",
         x = '',
         y = '') +
    theme(legend.position="top", legend.title = element_blank())
  
  gList[[2]] <- 
    ggplot(subset(data_All[data_All$fr == paste0('Failure rate: ', unique(fr)[2]),], t_ceiling <= 10), aes(x = t_ceiling, fill = group)) +  
    geom_histogram(aes(y = ..count..  / n), bins = 30) +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    geom_vline(xintercept=10, col = 'red', linetype="dashed") +
    facet_nested(fr ~ pi + z + cov, scales = 'free', nest_line = TRUE, bleed = TRUE) +
    labs(title = "",
         x = '',
         y = 'Proportion of Positive Tests') +
    theme(legend.position="none")
  
  gList[[3]] <- 
    ggplot(subset(data_All[data_All$fr == paste0('Failure rate: ', unique(fr)[3]),], t_ceiling <= 10), aes(x = t_ceiling, fill = group)) +  
    geom_histogram(aes(y = ..count..  / n), bins = 30) +
    scale_x_continuous(breaks = seq(0, 10, 1)) +
    geom_vline(xintercept=10, col = 'red', linetype="dashed") +
    facet_nested(fr ~ pi + z + cov, scales = 'free', nest_line = TRUE, bleed = TRUE) +
    labs(title = "",
         x = 'Age (days) at which test is administered',
         y = '') +
    theme(legend.position="none")  
  gAll <- marrangeGrob(gList, nrow=3, ncol=1, top = NULL)
  ggsave("scenarioII_V3.png", plot = gAll, width = 12, height = 10)
  
}

set.seed(20)
cov <- as.matrix(sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))) # binary variable, like gender, cov=0 as female, cov=1 as male
z <-  as.matrix(sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))) # Sample group indicator z with P(Z=0)=0.5

n <- 2500
pi0 <- head(rep(c(0.1, 0.2, 0.4), 3), -2)
failure_rate <- fr <- c(rep(0.9, 3), rep(0.5, 3), rep(0.3, 1))
value = matrix(round(c(-0.9, -4, 0.7, 0.1, 0.2, 2, 1.2,
                       -0.6, -3, 0.7, 0.1, 0.2, 2.3, 1.1,
                       2.2, -7, 0.7, 0.1, 0.3, 1.2, 0.9,
                       -0.9, -4, 0.7, 0.1, 0.2, 1.0, 0.1,
                       -0.1, -4, 0.7, 0.1, 0.2, 0.8, 0.1,
                       1.5, -4, 0.7, 0.3, 0.1, 0.2, 0.1,
                       -0.9, -4, 0.7, 0.1, 0.2, 0.1, 0.1), 2), byrow = TRUE, ncol = 7)
a0 <- value[,1]
aj <- value[,2]
gamma <- value[,c(4, 5)]
lambda <- value[,c(6, 7)]

simPlotsII(n, pi0, fr, a0, aj, beta, gamma, lambda, z, cov)

