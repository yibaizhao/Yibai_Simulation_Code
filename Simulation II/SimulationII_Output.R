source("SimulationII_Optim_EM.R")

################################ pi=0.1 failure_rate=0.9 #######################
pi0 = 0.1
failure_rate = 0.9

beta0 = 0.7
a0_0 = -0.9
aj_0 = -4
gamma0 = c(0.1, 0.2)
lambda0 = c(2, 1.2)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)

################################ pi=0.2 failure_rate=0.9 #######################
pi0 = 0.2
failure_rate = 0.9

beta0 = 0.7
a0_0 = -0.6
aj_0 = -3
gamma0 = c(0.1, 0.2)
lambda0 = c(2.3, 1.1)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)

################################ pi=0.4 failure_rate=0.9 #######################
pi0 = 0.4
failure_rate = 0.9

beta0 = 0.7
a0_0 = 2.2
aj_0 = -7
gamma0 = c(0.1, 0.3)
lambda0 = c(1.2, 0.9)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)

################################ pi=0.1 failure_rate=0.5 #######################
pi0 = 0.1
failure_rate = 0.5

beta0 = 0.7
a0_0 = -0.9
aj_0 = -4
gamma0 = c(0.1, 0.2)
lambda0 = c(1.0, 0.1)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)

################################ pi=0.2 failure_rate=0.5 #######################
pi0 = 0.2
failure_rate = 0.5

beta0 = 0.7
a0_0 = -0.1
aj_0 = -4
gamma0 = c(0.1, 0.2)
lambda0 = c(0.8, 0.1)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)

################################ pi=0.4 failure_rate=0.5 #######################
pi0 = 0.4
failure_rate = 0.5

beta0 = 0.7
a0_0 = 1.5
aj_0 = -4
gamma0 = c(0.3, 0.1)
lambda0 = c(0.2, 0.1)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)

################################ pi=0.1 failure_rate=0.3 #######################
pi0 = 0.1
failure_rate = 0.3

beta0 = 0.7
a0_0 = -0.9
aj_0 = -4
gamma0 = c(0.1, 0.2)
lambda0 = c(0.1, 0.1)


n = 2500

outF(n, failure_rate, pi0, a0_0, aj_0, beta0, gamma0, lambda0)
