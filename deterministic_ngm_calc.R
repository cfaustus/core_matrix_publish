############################################################
# Code for NGM R0
## FD & DD transmission
############################################################
# R0 at different phi + psi
# just by calculating a pathogen that comes into a completely naive population 
# at carrying capacity for that proportion forested
############################################################

rm(list = ls())

library(matrixcalc)
source('deterministic/ngm_cm_func.R')

# Range of phi
phi = seq(0.01, 1, by = 0.001) #vector for different forested proportions
epsilon = (1+cos(phi*(pi*3/2)-2.3))

btw= seq(0,1, by=0.05) # aka psi

############################################################
## DENSITY DEPENDENT R0

# Parameters defined (core; matrix) for DD transmission
paramsD <- list(d = c(0.1, 0.02),
                k = c(100, 100),
                beta.wDD = c(0.00136, 0.0004),
                gamma = c(0.03, 0.05),
                alpha = c(0.001, 0.01),
                sigma = c(0.05, 0.05))

R0C = 0.5
R0M = 1.5
(paramsD[['alpha']][1] + paramsD[['gamma']][1] + paramsD[['d']][1])*R0C/paramsD[['k']][1]
(paramsD[['alpha']][2] + paramsD[['gamma']][2] + paramsD[['d']][2])*R0M/paramsD[['k']][2]

R0_c = (paramsD[['beta.wDD']][1]*paramsD[['k']][1])/(paramsD[['alpha']][1] + paramsD[['gamma']][1] + paramsD[['d']][1])
R0_m = (paramsD[['beta.wDD']][2]*paramsD[['k']][2])/(paramsD[['alpha']][2] + paramsD[['gamma']][2] + paramsD[['d']][2])
print(c(R0_c,R0_m))

R0.c.dd = numeric(length(phi)) #assuming no matrix species present (or not transmissible)
R0.m.dd = numeric(length(phi)) #assuming no patch species present (or not transmissible)
R0.combin.dd = data.frame(matrix(NA, nrow = length(btw), ncol = length(phi)))

for (j in 1:length(btw)){
  for (i in 1:length(phi)){
    mat_F <- matrix.F.DD(paramsD[["beta.wDD"]], paramsD[["beta.wDD"]]*btw[j], epsilon[i], paramsD[['k']], phi[i])
    mat_V <- matrix.V(paramsD[["alpha"]], paramsD[["gamma"]], paramsD[["d"]])
    
    mat_G <- matrix.inverse(mat_V) %*% mat_F
    eigen.ngm.dd <- eigen(mat_G)
    R0.combin.dd[j,i] = eigen.ngm.dd$values[[1]]
    R0.c.dd[i] = (paramsD[["beta.wDD"]][1]*paramsD[["k"]][1]*(1.01-phi[i]))/
      (paramsD[["alpha"]][1] + paramsD[["gamma"]][1] + paramsD[["d"]][1])
    R0.m.dd[i] = (paramsD[["beta.wDD"]][2]*paramsD[["k"]][2]*phi[i])/
      (paramsD[["alpha"]][2] + paramsD[["gamma"]][2] + paramsD[["d"]][2])
  }
}

plotR0(phi, R0.c.dd, R0.m.dd, R0.combin.dd)

#export for plotting
dim(R0.combin.dd)
export = rbind(R0.combin.dd,R0.c.dd,R0.m.dd)
dim(export)
write.csv(export,'deterministic/output/scenario3_DD.csv')

############################################################
## FREQUENCY DEPENDENT R0
beta_c = R0_c*(paramsD[['alpha']][1] + paramsD[['gamma']][1] + paramsD[['d']][1])
beta_m = R0_m*(paramsD[['alpha']][2] + paramsD[['gamma']][2] + paramsD[['d']][2])
print(c(beta_c, beta_m))
# R0; plotting within and between host R0
paramsF <- list(d = c(0.1, 0.02),
                k = c(100, 100),
                beta.wFD = c(beta_c, beta_m),
                gamma = c(0.03, 0.05),
                alpha = c(0.001, 0.01),
                sigma = c(0.05, 0.05))

R0.m.fd = numeric(length(phi)) #
R0.c.fd = numeric(length(phi)) #
R0.combin.fd = data.frame(matrix(NA, nrow = length(btw), ncol = length(phi)))

for (j in 1:length(btw)){
  for (i in 1:length(phi)){
    mat_F <- matrix.F.FD(paramsF[["beta.wFD"]], paramsF[["beta.wFD"]]*btw[j], 
                         epsilon[i])
    mat_V <- matrix.V(paramsF[["alpha"]], paramsF[["gamma"]], paramsF[["d"]])
    
    mat_G <- matrix.inverse(mat_V) %*% mat_F
    eigen.ngm.fd <- eigen(mat_G)
    R0.combin.fd[j,i] = eigen.ngm.fd$values[[1]]
    
    R0.c.fd[i] = paramsF[["beta.wFD"]][1]/
      (paramsF[["alpha"]][1] + paramsF[["gamma"]][1] + paramsF[["d"]][1])
    R0.m.fd[i] = paramsF[["beta.wFD"]][2]/
      (paramsF[["alpha"]][2] + paramsF[["gamma"]][2] + paramsF[["d"]][2])
  }
}

#export for plotting
dim(R0.combin.fd)
plotR0(phi, R0.c.fd, R0.m.fd, R0.combin.fd)
export = rbind(R0.combin.fd,R0.c.fd,R0.m.fd)
dim(export)
write.csv(export,'deterministic/output/scenario3_FD.csv')
