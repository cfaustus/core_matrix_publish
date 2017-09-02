######### 
# Function to run ODEs  
# with random generation of parameters
# raina & christina
########


# cluster version of function
super_ode <- function(x,times){
  library(deSolve) #library for ode function
  library(pracma)
  #x=param.list[[40]]
  core.matrix.model <- function(Time, State, Parameters) {
    with(as.list(c(State, Parameters)), {
      N.c = S.c+I.c+R.c
      N.m = S.m+I.m+R.m
      
      dS.c = N.c*(rmax.c+d.c)*(1-N.c/(k.c*(1.01-phi))) - ((S.c*beta.c*I.c)/N.c^kappa+(epsilon*beta.m*psi*I.m*S.c)/(N.c+epsilon*N.m)^kappa) - d.c*S.c 
      dI.c = (beta.c*S.c*I.c)/N.c^kappa+(epsilon*beta.m*psi*S.c*I.m)/(N.c+epsilon*N.m)^kappa - I.c*(alpha.c*d.c+d.c+gamma.c)
      dR.c = I.c*gamma.c-R.c*d.c
      
      dS.m = N.m*(rmax.m+d.m)*(1-N.m/(k.m*phi)) - ((beta.m*S.m*I.m)/N.m^kappa+(epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - d.m*S.m 
      dI.m = ((beta.m*S.m*I.m)/N.m^kappa+(epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - I.m*(alpha.m*d.c+d.m+gamma.m)
      dR.m = I.m*gamma.m-R.m*d.m
      
      return(list(c(S.c, I.c, R.c, S.m, I.m, R.m)))
    })
  }
  initial.values = c(S.c=(1.01-x[['phi']])*x[['k.c']],I.c=1,R.c=0,S.m=x[['phi']]*x[['k.m']],I.m=0,R.m=0) 
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




