######### 
# Function to run ODEs  with varying forest cover
# raina & christina
########

s1<- read.csv("lhc_nomatrixv_fd.csv",header=T) # dataframe of random parameters
param.list <- apply(s1,1,as.list) #converting into a list
times <- seq(0, 50, by = 0.1) #time steps are constant for each run

# cluster version of function
super_ode <- function(x,times){
  library(deSolve) #library for ode function
  #x=param.list[[100]]
  patch.matrix.model <- function(Time, State, Parameters) {
    with(as.list(c(State, Parameters)), {
      N.p<-S.p+I.p+R.p
      N.m<-S.m+I.m+R.m
      
      dS.p<-N.p*(rmax.p+d.p)*(1-N.p/(k.p*(1.01-phi))) - ((S.p*beta.pp*I.p)/N.p^kappa+(epsilon*beta.mp*beta.pp*I.m*S.p)/(N.p+epsilon*N.m)^kappa) - d.p*S.p 
      dI.p<-(beta.pp*S.p*I.p)/N.p^kappa+(epsilon*beta.mp*beta.pp*S.p*I.m)/(N.p+epsilon*N.m)^kappa - I.p*(alpha.p*d.p+d.p+gamma.p)
      dR.p<-I.p*gamma.p-R.p*d.p
      
      dS.m<-N.m*(rmax.m+d.m)*(1-N.m/(k.m*phi)) - ((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.mm*beta.pm*S.m*I.p)/(epsilon*N.p+N.m)^kappa) - d.m*S.m 
      dI.m<-((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.mm*beta.pm*S.m*I.p)/(epsilon*N.p+N.m)^kappa) - I.m*(alpha.m*d.p+d.m+gamma.m)
      dR.m<-I.m*gamma.m-R.m*d.m
      
      return(list(c(dS.p, dI.p, dR.p,dS.m, dI.m, dR.m)))
    })
  }
  initial.values<-c(S.p=(1.01-x[['phi']])*x[['k.p']],I.p=1,R.p=0,S.m=x[['phi']]*x[['k.m']],I.m=0,R.m=0) 
  out <- as.data.frame(ode(func=patch.matrix.model,y=initial.values,parms=x,times=times, method = 'ode45'))
  #return(out)
  max=min(c(nrow(out),length(times)))
  out.vec=c(out[max,2:7],(out[max,3]/(out[max,2]+out[max,3]+out[max,4])),(out[max,6]/(out[max,5]+out[max,6]+out[max,7])),max)
  return(out.vec)
}


#super_ode(param.list[[100]],times)




