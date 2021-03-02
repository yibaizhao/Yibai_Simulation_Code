library(plyr)
library(readr)
library(reshape2)
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
setwd('/Volumes/GoogleDrive/My Drive/Work/HIV/code')
######################################## Scenario I ######################################## 
## True value
fr <- c(rep(0.9, 3), rep(0.5, 3), rep(0.3, 2))
pi <- head(rep(c(0.1, 0.2, 0.4), 3), -1)
data_True <- data.frame(Parameters = rep(c('a0', 'aj', 'beta', 'gamma', 'lambda'), length(pi)),
                        value = round(c(-0.9, -4, 0.7, 0.290028, 0.910090,
                                   -0.1, -4, 0.7, 0.550054, 0.450044,
                                   1.6, -4.4, 0.7, 3.000200e-02, 7.100700e-01,
                                   -0.9, -4, 0.7, 0.690068, 0.070006,
                                   -0.1, -4, 0.7, 0.570056, 0.070006,
                                   1.7, -4.5, 0.7, 2.500240e-01 , 5.000400e-02,
                                   -0.9, -4, 0.7, 0.530052, 0.050004,
                                   -4.2, 4, 0.7, 3.300320e-01, 5.000400e-02), 2),
                        group = rep(paste0('failure rate: ', fr, ', pi: ', pi), each = 5),
                        fr = paste0('failure rate: ', rep(fr, each = 5)),
                        pi = paste0('pi: ', rep(pi, each = 5)))
data_True$fr <- factor(data_True$fr, levels = unique(data_True$fr))

## Load Gu's results
mydir1 <- './ObjIII/Gu'
myfiles1 = list.files(path=mydir1, pattern="*.csv", full.names=TRUE)
data0 = ldply(myfiles1, function(x){
  dt0 <- read_csv(x)[,-1]
  dt <- as.data.frame(t(matrix(as.numeric(dt0), ncol = 3)))
  dt$parameters <- c('beta', 'gamma', 'lambda')
  dt
  }
  )
colnames(data0) <- c('Estimates', 'Std', 'CovProb', 'Parameters')
data0$fr = paste0('failure rate: ', rep(c(rep(0.3, 2), rep(0.5, 3), rep(0.9, 3)), each = 3))
data0$pi = paste0('pi: ', rep(c(0.1, 0.2, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4), each = 3))
data <- merge(data0, data_True, by = c('Parameters', 'fr', 'pi'))
data$`Relative Bias` <- (data$Estimates - data$value) / data$value
data$group <- rep(paste0('failure rate: ', 
                         c(rep(0.3, 2), rep(0.5, 3), rep(0.9, 3)), 
                         ', pi: ',
                         c(0.1, 0.2, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(data) / length(myfiles1))
# data$fr <- paste0('failure rate: ', rep(c(rep(0.3, 2), rep(0.5, 3), rep(0.9, 3)), each = 3))
# data$pi <- paste0('pi: ', rep(c(0.1, 0.2, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4), each = 3))
data_Gu <- melt(data, id.vars = colnames(data)[c(1:3, 8)])
data_Gu$method <- "Non-Mixture Model"
## Load proposed outputs
mydir2 <- './ObjIII/Numerical_Method/Results'
myfiles2 = list.files(path=mydir2, pattern="*.csv", full.names=TRUE)
data0 = ldply(myfiles2, read_csv)
names(data0)[1] <- names(data_True)[1]
data0$fr = paste0('failure rate: ', rep(c(rep(0.3, 2), rep(0.5, 3), rep(0.9, 3)), each = 5))
data0$pi = paste0('pi: ', rep(c(0.1, 0.2, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4), each = 5))
data0$Parameters[data0$Parameters == 'a1'] <- 'aj'
data <- merge(data0, data_True, by = c('Parameters', 'fr', 'pi'))
data$`Relative Bias` <- (data$Estimates - data$value) / data$value
data$group <- rep(paste0('failure rate: ', 
                     c(rep(0.3, 2), rep(0.5, 3), rep(0.9, 3)), 
                     ', pi: ',
                     c(0.1, 0.2, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(data) / length(myfiles2))
data_proposed <- melt(data, id.vars = colnames(data)[c(1, 2:3, 8)])
data_proposed$method <- 'Mixture Model'
data_proposed$Parameters[data_proposed$Parameters == 'a1'] <- 'aj'

data <- do.call(rbind, list(data_Gu, data_proposed))
data$fr <- factor(data$fr, levels = rev(unique(data$fr)))

library(ggplot2)
data_Est <- data[data$variable == 'Estimates',]
ggplot(data_Est, aes(x = Parameters, y = value, color = Parameters)) +
  geom_point(aes(shape = method, group = group), alpha = 0.5, size=5, position=position_jitter(w=0.15)) + #, position=position_jitter(h=0.15,w=0.15))
  geom_point(data = data_True, shape = '*', size = 10, alpha = 0.8, position=position_jitter(w=0.15)) +
  labs(y = "Estimates", caption = "* refers to the true value") +
  facet_grid(fr ~ pi) +
  scale_x_discrete(labels = expression(alpha[0], alpha[1], beta, gamma, lambda)) +
  scale_color_discrete(labels = expression(alpha[0], alpha[1], beta, gamma, lambda)) +
  theme_bw()
ggsave('ScenarioI_Est_V4.png', width = 8, height = 12)

data_relative_bias <- data[data$variable == 'Relative Bias',]
ggplot(data_relative_bias, aes(x = Parameters, y = value, shape = method, group = group, color = Parameters)) +
  geom_point(alpha = 0.8, size=5, position=position_jitter(w=0.15))+
  labs(y = "Relative Bias") +
  facet_grid(fr ~ pi) +
  scale_x_discrete(labels=expression(alpha[0], alpha[1], beta, gamma, lambda)) +
  scale_color_discrete(labels = expression(alpha[0], alpha[1], beta, gamma, lambda)) +
  theme_bw()
ggsave('ScenarioI_RB_V4_Slides.png', width = 10, height = 7)


data_CovProb <- data[data$variable == 'CovProb',]
ggplot(data_CovProb, aes(x = Parameters, y = value, shape = method, group = group, color = Parameters)) +
  geom_point(alpha = 0.8, size=5, position=position_jitter(w=0.15)) +
  geom_hline(yintercept=0.95, col = 'red', linetype="dashed") + 
  scale_y_continuous(breaks = sort(c(seq(min(data_CovProb$value), max(data_CovProb$value), length.out=5), 0.95))) + 
  ylab('Coverage Probability') +
  facet_grid(fr ~ pi) +
  scale_x_discrete(labels=expression(alpha[0], alpha[1], beta, gamma, lambda)) +
  scale_color_discrete(labels = expression(alpha[0], alpha[1], beta, gamma, lambda)) +
  theme_bw()
ggsave('ScenarioI_CovProb_V4_Slides.png', width = 10, height = 7)


######################################## Scenario II ######################################## 
## True value
fr <- c(rep(0.9, 3), rep(0.5, 3), 0.3)
pi <- head(rep(c(0.1, 0.2, 0.4), 3), -2)
data_True <- data.frame(Parameters = rep(c('a0', 'aj', 'beta', 'gamma1', 'gamma2', 'lambda1', 'lambda2'), length(pi)),
                        value = round(c(-0.9, -4, 0.7, 0.1, 0.2, 2, 1.2,
                                -0.6, -3, 0.7, 0.1, 0.2, 2.3, 1.1,
                                2.2, -7, 0.7, 0.1, 0.3, 1.2, 0.9,
                                -0.9, -4, 0.7, 0.1, 0.2, 1.0, 0.1,
                                -0.1, -4, 0.7, 0.1, 0.2, 0.8, 0.1,
                                1.5, -4, 0.7, 0.3, 0.1, 0.2, 0.1,
                                -0.9, -4, 0.7, 0.1, 0.2, 0.1, 0.1), 2),
                        group = rep(paste0('failure rate: ', fr, ', pi: ', pi), each = 7),
                        fr = paste0('failure rate: ', rep(fr, each = 7)),
                        pi = paste0('pi: ', rep(pi, each = 7)))
data_True$fr <- factor(data_True$fr, levels = unique(data_True$fr))


## Load Gu's results
mydir1 <- './ObjIV/Gu'
myfiles1 = list.files(path=mydir1, pattern="*.csv", full.names=TRUE)
data0 = ldply(myfiles1, function(x){
  dt0 <- as.data.frame(read_csv(x))[,2]
  dt <- as.data.frame((matrix(dt0, byrow = TRUE, ncol = 3)))
  dt$parameters <- c('beta', 'gamma1', 'gamma2', 'lambda1', 'lambda2')
  dt
}
)
colnames(data0) <- c('Estimates', 'Std', 'CovProb', 'Parameters')
data0$fr = paste0('failure rate: ', rep(c(0.3, rep(0.5, 3), rep(0.9, 3)), each = 5))
data0$pi = paste0('pi: ', rep(c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4), each = 5))
data <- merge(data0, data_True, by = c('Parameters', 'fr', 'pi'))
data$`Relative Bias` <- (data$Estimates - data$value) / data$value

data$group <- rep(paste0('failure rate: ', 
                         c(0.3, rep(0.5, 3), rep(0.9, 3)), 
                         ', pi: ',
                         c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(data) / length(myfiles1))
data_Gu <- melt(data, id.vars = colnames(data)[c(1:3, 8)])
data_Gu$Method <- "Stratified Non-Mixture Model"
## Load proposed outputs
mydir2 <- './ObjIV/Numeric_Method/Results'
myfiles2 = list.files(path=mydir2, pattern="*.csv", full.names=TRUE)
data0 = ldply(myfiles2, read_csv)
colnames(data0)[1] <- 'Parameters'
data0$fr = paste0('failure rate: ', rep(c(0.3, rep(0.5, 3), rep(0.9, 3)), each = 7))
data0$pi = paste0('pi: ', rep(c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4), each = 7))
data0$Parameters[data0$Parameters == 'a1'] <- 'aj'
data <- merge(data0, data_True, by = c('Parameters', 'fr', 'pi'))
data$`Relative Bias` <- (data$Estimates - data$value) / data$value

data$group <- rep(paste0('failure rate: ', 
                         c(0.3, rep(0.5, 3), rep(0.9, 3)),
                         ', pi: ',
                         c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(data) / length(myfiles2))
data_proposed <- melt(data, id.vars = colnames(data)[c(1:3, 8)])
data_proposed$Method <- 'Stratified Mixture Model (Proposed)'

data <- do.call(rbind, list(data_Gu, data_proposed))
data$fr <- factor(data$fr, levels = rev(unique(data$fr)))

library(ggplot2)
data_Est <- data[data$variable == 'Estimates',]
ggplot(data_Est, aes(x = Parameters, y = value, color = Parameters)) +
  geom_point(aes(shape = Method, group = group), alpha = 0.5, size=5, position=position_jitter(w=0.15)) + #, position=position_jitter(h=0.15,w=0.15))
  geom_point(data = data_True, shape = '*', size = 10, alpha = 0.8, position=position_jitter(w=0.15)) +
  labs(y = "Estimates", caption = "* refers to the true value") +
  facet_grid(fr ~ pi) +
  scale_x_discrete(labels=expression(alpha[0], alpha[1], beta, gamma[1], gamma[2], lambda[1], lambda[2])) +
  scale_color_discrete(labels = expression(alpha[0], alpha[1], beta, gamma[1], gamma[2], lambda[1], lambda[2])) +
  theme_bw()
ggsave('ScenarioII_Est_V4.png', width = 8, height = 12)

data_relative_bias <- data[data$variable == 'Relative Bias',]
ggplot(data_relative_bias, aes(x = Parameters, y = value, shape = Method, group = group, color = Parameters)) +
  geom_point(alpha = 0.8, size=5, position=position_jitter(w=0.15))+
  labs(y = "Relative Bias") +
  facet_grid(fr ~ pi) +
  scale_x_discrete(labels=expression(alpha[0], alpha[1], beta, gamma[1], gamma[2], lambda[1], lambda[2])) +
  scale_color_discrete(labels = expression(alpha[0], alpha[1], beta, gamma[1], gamma[2], lambda[1], lambda[2])) +
  theme_bw()
ggsave('ScenarioII_RB_V4_Slides.png', width = 10, height = 7)


data_CovProb <- data[data$variable == 'CovProb',]
ggplot(data_CovProb, aes(x = Parameters, y = value, shape = Method, group = group, color = Parameters)) +
  geom_point(alpha = 0.8, size=5, position=position_jitter(w=0.15)) +
  geom_hline(yintercept=0.95, col = 'red', linetype="dashed") + 
  scale_y_continuous(breaks = sort(c(seq(min(data_CovProb$value), max(data_CovProb$value), length.out=5), 0.95))) + 
  ylab('Coverage Probability') +
  facet_grid(fr ~ pi) +
  scale_x_discrete(labels=expression(alpha[0], alpha[1], beta, gamma[1], gamma[2], lambda[1], lambda[2])) +
  scale_color_discrete(labels = expression(alpha[0], alpha[1], beta, gamma[1], gamma[2], lambda[1], lambda[2])) +
  theme_bw()
ggsave('ScenarioII_CovProb_V4_Slides.png', width = 10, height = 7)

#============================= Generate Tables of Scenario I ============================= 
## True value
fr <- c(rep(0.9, 3), rep(0.5, 3), rep(0.3, 2))
pi <- head(rep(c(0.1, 0.2, 0.4), 3), -1)
dataI_True <- data.frame(Parameters = rep(c('a0', 'aj', 'beta', 'gamma', 'lambda'), length(pi)),
                        value = round(c(-0.9, -4, 0.7, 0.290028, 0.910090,
                                        -0.1, -4, 0.7, 0.550054, 0.450044,
                                        1.6, -4.4, 0.7, 3.000200e-02, 7.100700e-01,
                                        -0.9, -4, 0.7, 0.690068, 0.070006,
                                        -0.1, -4, 0.7, 0.570056, 0.070006,
                                        1.7, -4.5, 0.7, 2.500240e-01 , 5.000400e-02,
                                        -0.9, -4, 0.7, 0.530052, 0.050004,
                                        -4.2, 4, 0.7, 3.300320e-01, 5.000400e-02), 2),
                        group = rep(paste0('failure rate: ', fr, ', pi: ', pi), each = 5),
                        fr = paste0('failure rate: ', rep(fr, each = 5)),
                        pi = paste0('pi: ', rep(pi, each = 5)))
dataI_True$fr <- factor(dataI_True$fr, levels = unique(dataI_True$fr))
dataI_True$from <- 'True'
dataI_Truel <- split(dataI_True, dataI_True$Parameters)
## Load proposed outputs
mydir2 <- './ObjIII/Numerical_Method/Results'
myfiles2 = list.files(path=mydir2, pattern="*.csv", full.names=TRUE)
dataI_Est = ldply(myfiles2, read_csv)
dataI_Est$group <- rep(paste0('failure rate: ', 
                         c(rep(0.3, 2), rep(0.5, 3), rep(0.9, 3)), 
                         ', pi: ',
                         c(0.1, 0.2, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(dataI_Est) / length(myfiles2))
colnames(dataI_Est)[1] <- 'Parameters'
dataI_Est$group <- factor(dataI_Est$group, levels = unique(dataI_True$group))
dataI_Est <- dataI_Est[order(dataI_Est$group),]
dataI_Est$from <- 'Ests'
dataI_Est$Parameters[dataI_Est$Parameters == 'a1'] <- 'aj'
dataI_Estl <- split(dataI_Est, unique(dataI_Est$Parameters))

dataIl <- lapply(dataI_Estl, function(x){
  # Data True
  data_True <- dataI_Truel[[unique(x$Parameters)]]
  true_val <- round(data_True$value, 2)
  # Data Ests
  std <- paste0(round(x$Estimates, 2), ' (', round(x$Std, 2), ')')
  dt <- c(unique(x$Parameters), std)
  dt2 <- as.data.frame(rbind(dt, c(paste0('True ', unique(data_True$Parameters)), true_val), c('Cov.Prob', round(x$CovProb, 2))))     
  colnames(dt2) <- c('Parameters', as.character(x$group))
  dt2
})
dataI <- do.call(rbind, dataIl)
write.csv(dataI, 'dataI_V2.csv', row.names = FALSE)


#============================= Generate Tables of Scenario II ============================= 
## True value
fr <- c(rep(0.9, 3), rep(0.5, 3), 0.3)
pi <- head(rep(c(0.1, 0.2, 0.4), 3), -2)
dataII_True <- data.frame(Parameters = rep(c('a0', 'aj', 'beta', 'gamma1', 'gamma2', 'lambda1', 'lambda2'), length(pi)),
                        value = round(c(-0.9, -4, 0.7, 0.1, 0.2, 2, 1.2,
                                        -0.6, -3, 0.7, 0.1, 0.2, 2.3, 1.1,
                                        2.2, -7, 0.7, 0.1, 0.3, 1.2, 0.9,
                                        -0.9, -4, 0.7, 0.1, 0.2, 1.0, 0.1,
                                        -0.1, -4, 0.7, 0.1, 0.2, 0.8, 0.1,
                                        1.5, -4, 0.7, 0.3, 0.1, 0.2, 0.1,
                                        -0.9, -4, 0.7, 0.1, 0.2, 0.1, 0.1), 2),
                        group = rep(paste0('failure rate: ', fr, ', pi: ', pi), each = 7),
                        fr = paste0('failure rate: ', rep(fr, each = 7)),
                        pi = paste0('pi: ', rep(pi, each = 7)))
dataII_True$fr <- factor(dataII_True$fr, levels = unique(dataII_True$fr))
dataII_Truel <- split(dataII_True, dataII_True$Parameters)


## Load Gu's results
mydir1 <- './ObjIV/Gu'
myfiles1 = list.files(path=mydir1, pattern="*.csv", full.names=TRUE)
data = ldply(myfiles1, function(x){
  dt0 <- as.data.frame(read_csv(x))[,2]
  dt <- as.data.frame((matrix(dt0, byrow = TRUE, ncol = 3)))
  dt$parameters <- c('beta', 'gamma1', 'gamma2', 'lambda1', 'lambda2')
  dt
}
)
colnames(data) <- c('Estimates', 'Std', 'CovProb', 'Parameters')
data$group <- rep(paste0('failure rate: ', 
                         c(0.3, rep(0.5, 3), rep(0.9, 3)), 
                         ', pi: ',
                         c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(data) / length(myfiles1))
data$fr <- paste0('failure rate: ', rep(c(0.3, rep(0.5, 3), rep(0.9, 3)), each =  nrow(data) / length(myfiles1)))
data$pi <- paste0('pi: ', rep(c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4), each =  nrow(data) / length(myfiles1)))
data_Gu <- melt(data, id.vars = colnames(data)[c(4, 5, 6, 7)])
data_Gu$method <- "Gu's Stratified Non-Mixture Model"
## Load proposed outputs
mydir2 <- './ObjIV/Numeric_Method/Results'
myfiles2 = list.files(path=mydir2, pattern="*.csv", full.names=TRUE)
dataII_Est = ldply(myfiles2, read_csv)
dataII_Est$group <- rep(paste0('failure rate: ', 
                         c(0.3, rep(0.5, 3), rep(0.9, 3)),
                         ', pi: ',
                         c(0.1, 0.1, 0.2, 0.4, 0.1, 0.2, 0.4)),
                  each = nrow(dataII_Est) / length(myfiles2))
colnames(dataII_Est)[1] <- 'Parameters'
dataII_Est$group <- factor(dataII_Est$group, levels = unique(dataII_True$group))
dataII_Est <- dataII_Est[order(dataII_Est$group),]
dataII_Estl <- split(dataII_Est, unique(dataII_Est$Parameters))

dataIIl <- lapply(dataII_Estl, function(x){
  # Data True
  data_True <- dataII_Truel[[unique(x$Parameters)]]
  true_val <- round(data_True$value, 2)
  # Data Ests
  std <- paste0(round(x$Estimates, 2), ' (', round(x$Std, 2), ')')
  dt <- c(unique(x$Parameters), std)
  dt2 <- as.data.frame(rbind(dt, c(paste0('True ', unique(data_True$Parameters)), true_val), c('Cov.Prob', round(x$CovProb, 2))))     
  colnames(dt2) <- c('Parameters', as.character(x$group))
  dt2
})
dataII <- do.call(rbind, dataIIl)
write.csv(dataII, 'dataII_V2.csv', row.names = FALSE)

