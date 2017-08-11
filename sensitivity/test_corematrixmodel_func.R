core.matrix.model <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.c<-S.c+I.c+R.c
    N.m<-S.m+I.m+R.m
    
    dS.c<-N.c*(rmax.c+d.c)*(1-N.c/(k.c*(1.01-phi))) - ((S.c*beta.c*I.c)/N.c^kappa+(epsilon*psi*beta.m*I.m*S.c)/(N.c+epsilon*N.m)^kappa) - d.c*S.c 
    dI.c<-(beta.c*S.c*I.c)/N.c^kappa+(epsilon*psi*beta.m*S.c*I.m)/(N.c+epsilon*N.m)^kappa - I.c*(alpha.c*d.c+d.c+gamma.c)
    dR.c<-I.c*gamma.c-R.c*d.c
    
    dS.m<-N.m*(rmax.m+d.m)*(1-N.m/(k.m*phi)) - ((beta.m*S.m*I.m)/N.m^kappa+(epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - d.m*S.m 
    dI.m<-((beta.m*S.m*I.m)/N.m^kappa+(epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - I.m*(alpha.m*d.c+d.m+gamma.m)
    dR.m<-I.m*gamma.m-R.m*d.m
    
    return(list(c(dS.c, dI.c, dR.c,dS.m, dI.m, dR.m)))
  })
}

core.matrix.model_f <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.c<-S.c+I.c+R.c
    N.m<-S.m+I.m+R.m
    
    dS.c<-N.c*(rmax.c+d.c)*(1-N.c*k.c/(1.01-phi)) - ((S.c*beta.cp*I.c)/N.c^kappa+(epsilon*psi*beta.m*I.m*S.c)/(N.c+epsilon*N.m)^kappa) - d.c*S.c 
    dI.c<-(beta.c*S.c*I.c)/N.c^kappa+(epsilon*psi*beta.m*S.c*I.m)/(N.c+epsilon*N.m)^kappa - I.c*(alpha.c*d.c+d.c+gamma.c)
    dR.c<-I.c*gamma.c-R.c*d.c
    
    dS.m<-N.m*(rmax.m+d.m)*(1-N.m*k.m/phi) - ((beta.m*S.m*I.m)/N.m^kappa+(epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - d.m*S.m 
    dI.m<-((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*psi*beta.c*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - I.m*(alpha.m*d.c+d.m+gamma.m)
    dR.m<-I.m*gamma.m-R.m*d.m
    
    return(list(c(dS.c, dI.c, dR.c,dS.m, dI.m, dR.m)))
  })
}