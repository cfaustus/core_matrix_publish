###########################################################3
###
###     Aaron - Gillespie Algorithm for SIR epidemic
###     Modified by Hamish to add spillover
###     Density dependent transmission
###     Parameterized for Ebola Uganda 2000-2001
###     Case mortality 53%
###     mean time to death 8 days => death rate/year 365/8= 42.625
###     case mortality rate 53%   => recovery rate          40.46
###
###
###     R0=1.79 (from Ferrari, M.J., Bjornstad, O.N. & Dobson, A.P. (2005) Mathematical Biosciences, 198, 14-26.
###     beta= R0*(alpha+mu+gamma)
###     beta=148.7
###  First set routine to define one-step function

SIR.onestep <- function (x, params) {
  X <- x[2]
  Y <- x[3]
  Z <- x[4]
  N <- X+Y+Z
  beta <- params["beta"]
  mu <- params["mu"]
  gamma <- params["gamma"]
  spill<-params["spill"]
  alpha<-params["alpha"]
  
  ## each individual rate
  
  rates <- c(
    birth=mu*N,
    infection=beta*X*Y/N+spill*X,
    recovery=gamma*Y,
    sdeath=mu*X,
    ideath=(mu+alpha)*Y,
    rdeath=mu*Z
  )
  
  #### what changes with each event?
  
  transitions <- list( 
    birth=c(1,0,0),
    infection=c(-1,1,0),
    recovery=c(0,-1,1),
    sdeath=c(-1,0,0),
    ideath=c(0,-1,0),
    rdeath=c(0,0,-1)
  )
  
  ## total event rate
  total.rate <- sum(rates)
  
  ## waiting time (note exponentially distributed random events)
  
  tau <- rexp(n=1,rate=total.rate)
  
  ## which event occurs?
  event <- sample.int(n=6,size=1,prob=rates/total.rate)
  
  x+c(tau,transitions[[event]],as.numeric(event==2))## double square bracket as recursive from list
#number of cases incremented by 1 if event 2 (infection) occurs  
}

############   Routine to define bounds of simulation

SIR.simul <- function (x, params, maxstep = 10000) {
  output <- array(dim=c(maxstep+1,5))
  colnames(output) <- names(x)
  output[1,] <- x
  k <- 1
  ## loop until either k > maxstep or
  ## time >1 (note that the loop will COMMENCE unless t>1, so final time >1 )
  while ((k <= maxstep) && (x[1] < 1)) {
    k <- k+1
    output[k,] <- x <- SIR.onestep(x,params)
    
  }
  as.data.frame(output[1:k,])
}

#####   Routine to run program  ~ note use of plyr as rdply to store results
require(plyr) ###  plyr really powerful for simplifying repeated complex operations
set.seed(56856583)
spillover<-seq(0.002,0.04,by=0.002)
breaks<-seq(0,400,length.out=201) #sets up bins by 2 up to 200
out.matrix<-matrix(NA,nrow=length(spillover),ncol=length(breaks)-1)

for(i in 1:length(spillover)){
  nsims <- 1000
  xstart <- c(time=0,X=400,Y=0,Z=0,cases=0) #initial conditions
  params <- c(mu=0.02,beta=148.7,gamma=40.46,spill=spillover[i],alpha=42.62) #parameters
  
  
  simdat <- rdply(
    nsims,
    SIR.simul(xstart,params)
  )
  results<-aggregate(cases~.n,data=simdat,FUN=max) #finds maximum number of cases per iteration
  out.hist<-hist(results$cases,plot=F,breaks=breaks) #doesn't plot, but groups cases into bins
  out.matrix[i,]<-out.hist$counts #puts counts into out.matrix
  rm(simdat,results,out.hist) #clean up after each 1000 iterations 
  
}
library(RColorBrewer)
image(z=log10(out.matrix+1),x=spillover,y=breaks,xlab="spillover",ylab="cases",col=brewer.pal(8,"YlOrRd"))
pdf("ebola freq dep.pdf")
image(z=log10(out.matrix+1),x=spillover,y=breaks,xlab="spillover",ylab="cases",col=brewer.pal(8,"YlOrRd"))
dev.off()



