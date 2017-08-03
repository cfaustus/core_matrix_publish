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

patch.matrix.model_f <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.p<-S.p+I.p+R.p
    N.m<-S.m+I.m+R.m
    
    dS.p<-N.p*(rmax.p+d.p)*(1-N.p*k.p/f) - ((S.p*beta.pp*I.p)/N.p^kappa+(epsilon*beta.mp*beta.pp*I.m*S.p)/(N.p+epsilon*N.m)^kappa) - d.p*S.p 
    dI.p<-(beta.pp*S.p*I.p)/N.p^kappa+(epsilon*beta.mp*beta.pp*S.p*I.m)/(N.p+epsilon*N.m)^kappa - I.p*(alpha.p*d.p+d.p+gamma.p)
    dR.p<-I.p*gamma.p-R.p*d.p
    
    dS.m<-N.m*(rmax.m+d.m)*(1-N.m*k.m/(1.5-f)) - ((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.mm*beta.pm*S.m*I.p)/(epsilon*N.p+N.m)^kappa) - d.m*S.m 
    dI.m<-((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.mm*beta.pm*S.m*I.p)/(epsilon*N.p+N.m)^kappa) - I.m*(alpha.m*d.p+d.m+gamma.m)
    dR.m<-I.m*gamma.m-R.m*d.m
    
    return(list(c(dS.p, dI.p, dR.p,dS.m, dI.m, dR.m)))
  })
}