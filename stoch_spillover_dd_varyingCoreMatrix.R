###########################################################
###     Aaron - Gillespie Algorithm for SIR epidemic
###     Modified by Hamish to add spillover
###     modified to include core-matrix variations by C Faust
###     Density dependent transmission

rm(list = ls())

#libraries
library(plyr) 
library(RColorBrewer)

N_sim = 400 #arbitrary
R0_sim = 1.79 # R0 from Ferrari et al. 2005
mu_sim = 0.02 # est mean lifespan of 50 years
alpha_sim = 45.62 # mean time to death 8 days (365/8) CDC 2001
gamma_sim = 40.46 # case mortality rate 53% 
beta_sim = R0_sim*(mu_sim + gamma_sim + alpha_sim)/N_sim # 0.3852975

source('spillover_dd_func.R')
# functions: SIR.simul and SIR.onestep

#####   Routine to run program 
set.seed(56856583)

### how parameters change over forestation gradient
phi = seq(0,1, by=0.01)
K = 1+phi*400
K.c = 400*(1.01-phi)
epsilon = (1+cos(phi*(pi*3/2)-2.3)) 
spillover = K.c * 0.0001 *epsilon
plot(phi, spillover, pch=20, ylab= "spillover force of infection",
     main = "with edge effects", bty = 'n')
plot(phi, epsilon, ylim=c(0,2))
lines(phi, spillover, col='darkred')
mat <- as.data.frame(matrix(NA, nrow=length(phi), ncol=5))
epidemic.size<-seq(0, 400, length.out = 401) #sets up bins by 2 up to 200
nsims=1000
out.matrix<-matrix(NA,nrow=length(phi),ncol=length(epidemic.size)-1)
cases.matrix<-matrix(NA,nrow=length(phi),ncol=nsims)

for(i in 2:length(phi)){
  #i=2
  #initial conditions
  xstart <- c(time=0,
              X=K[i],
              Y=0,
              Z=0,
              cases=0) 
  params <- c(mu = mu_sim,
              beta = beta_sim,
              gamma = gamma_sim,
              spill=spillover[i],
              alpha = alpha_sim) #parameters
  #run silumations
  simdat <- rdply(
    nsims,
    SIR.simul(xstart,params)
  )
  results<-aggregate(cases~.n,data=simdat,FUN=max) #finds maximum number of cases per iteration
  cases.matrix[i,]<-results$cases
  out.hist<-hist(results$cases,plot=T,breaks=epidemic.size) #doesn't plot, but groups cases into bins
  out.matrix[i,]<-out.hist$counts #puts counts into out.matrix
  rm(simdat,results,out.hist) #clean up after each 1000 iterations 
}

dim(out.matrix)
summary(out.matrix)
out<-as.data.frame(out.matrix)
colnames(out) <-epidemic.size[0:400]
rownames(out)<- phi
write.csv(out, file = "output/ddsimulationwithepsilon_phi_400sims_24may.csv", row.names=TRUE)

dim(cases.matrix)
summary(cases.matrix)
out2<-as.data.frame(cases.matrix)
rownames(out2)<- phi
write.csv(out2, file = "output/ddsimulationwithepsilon_phi_400sims_cases_24may.csv", row.names=TRUE)
