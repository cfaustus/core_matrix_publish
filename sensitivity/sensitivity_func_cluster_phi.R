######### 
# Function to run ODEs  
# with random generation of parameters
# raina & christina
########
s1<- read.csv("sensitivity/lhc/lhc_humans_fd_R0vary_phi0.1.csv",header=T) # dataframe of random parameters
param.list <- apply(s1,1,as.list) #converting into a list
times <- seq(0, 20, by = 0.01) #time steps are constant for each run

# cluster version of function
super_ode <- function(x,times){
  library(deSolve) #library for ode function
  library(pracma)
  #x=param.list[[1000]]
  core.matrix.model <- function(Time, State, Parameters) {
    with(as.list(c(State, Parameters)), {
      N.c = S.c+I.c+R.c
      N.m = S.m+I.m+R.m

      dS.c <- ifelse((1 - N.c/(k.c*(1.01-phi))) > 0, #IF UNDER CARRYING CAPACITY THEN..
                     ((N.c)*(rmax.c+d.c)*(1 - N.c/(k.c*(1.01-phi))) - ((beta.c*S.c*I.c)/N.c^kappa
                        + (epsilon*psi*beta.m*S.c*I.m)/(N.c+epsilon*N.m)^kappa)
                        - d.c*S.c + sigma.c*R.c),
                     # ELSE... NO BIRTHS
                     (- ((beta.c*S.c*I.c)/N.c^kappa + (epsilon*psi*beta.m*S.c*I.m)/(N.c+epsilon*N.m)^kappa) 
                      - d.c*S.c + sigma.c*R.c)) 
      dI.c <- (beta.c*S.c*I.c)/N.c^kappa + (epsilon*phi*beta.m*S.c*I.m)/(N.c+epsilon*N.m)^kappa  - I.c*(alpha.c*d.c+d.c+gamma.c)
      dR.c <- I.c*gamma.c - R.c*(d.c+sigma.c)
      
      dS.m <-  ifelse((1 - N.m/(k.m*phi)) > 0,
                      ((N.m)*(rmax.m+d.m)*(1 - N.m/(k.m*phi)) - ((beta.m*S.m*I.m)/N.m^kappa
                         + (epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) 
                         - d.m*S.m  + sigma.m*R.m), 
                      - ((beta.m*S.m*I.m)/N.m^kappa + (epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - d.m*S.m + sigma.m*R.m)
      dI.m <- ((beta.m*S.m*I.m)/N.m^kappa + (epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - I.m*(alpha.m*d.m+d.m+gamma.m)
      dR.m <- I.m*gamma.m - R.m*(d.m+sigma.m)
      
      return(list(c(dS.c, dI.c, dR.c,dS.m, dI.m, dR.m)))
    })
  }
  initial.values = c(S.c=(1.01-x[['phi']])*x[['k.c']]*0.9,I.c=1,R.c=0,S.m=x[['phi']]*x[['k.m']]*0.9,I.m=0,R.m=0) 
  out = as.data.frame(ode(func=core.matrix.model,y=initial.values,parms=x,times=times, method = 'ode45'))
  #return(out)
  max=min(c(nrow(out),length(times)))
  out.vec=c(out[max,2:7],
            (out[max,3]/(out[max,2]+out[max,3]+out[max,4])),
            (out[max,6]/(out[max,5]+out[max,6]+out[max,7])),
            max,
            trapz(out$time,out[,3]/(out[,2]+out[,3]+out[,4])),
            trapz(out$time,out[,6]/(out[,5]+out[,6]+out[,7])),
            max(out[,3]/(out[,2]+out[,3]+out[,4])),
            max(out[,6]/(out[,5]+out[,6]+out[,7]))
            )
  return(out.vec)
}


#super_ode(param.list[[100]],times)




