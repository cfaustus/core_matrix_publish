patch.matrix.model.phi <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    
    N.c <- S.c+I.c+R.c
    N.m <- S.m+I.m+R.m
    #Dynamics in the Core Population
    dS.c <- ifelse((1 - N.c/(k.c*(1.01-phi))) > 0, #IF UNDER CARRYING CAPACITY THEN..
                   (N.c)*(b.c)*(1 - N.c/(k.c*(1.01-phi))) 
                   - ((beta.cc*S.c*I.c)/N.c^kappa
                      + (epsilon*beta.mc*beta.mm*S.c*I.m)/(N.c+epsilon*N.m)^kappa)
                   - d.c*S.c 
                   + sigma.c*R.c, # ELSE... NO BIRTHS
                   #- ((beta.cc*S.c*I.c)/N.c^kappa
                   #   + (epsilon*beta.mc*S.c*I.m)/(N.c+epsilon*N.m)^kappa)
                   - d.c*S.c 
                   + sigma.c*R.c) 
    
    dI.c <- (beta.cc*S.c*I.c)/N.c^kappa + (epsilon*beta.mc*beta.mm*S.c*I.m)/(N.c+epsilon*N.m)^kappa  - I.c*(alpha.c+d.c+gamma.c)
    dR.c <- I.c*gamma.c - R.c*(d.c+sigma.c)
    #Dynamics in the matrix
    dS.m <-  ifelse((1 - N.m/(k.m*phi)) > 0,
                    (N.m)*b.m*(1 - N.m/(k.m*phi)) 
                    - ((beta.mm*S.m*I.m)/N.m^kappa
                       + (epsilon*beta.cm*beta.cc*S.m*I.c)/(epsilon*N.c+N.m)^kappa) 
                    - d.m*S.m 
                    + sigma.m*R.m, 
                    #- ((beta.mm*S.m*I.m)/N.m^kappa
                    #    + (epsilon*beta.cm*S.m*I.c)/(epsilon*N.c+N.m)^kappa) 
                    - d.m*S.m
                    + sigma.m*R.m)
    dI.m <- ((beta.mm*S.m*I.m)/N.m^kappa + (epsilon*beta.cm*beta.cc*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - I.m*(alpha.m+d.m+gamma.m)
    dR.m <- I.m*gamma.m - R.m*(d.m+sigma.m)
    
    return(list(c(dS.c, dI.c, dR.c,dS.m, dI.m, dR.m)))
  })
}

patch.matrix.model.noifthen <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.c<-S.c+I.c+R.c
    N.m<-S.m+I.m+R.m
    
    dS.c<-(N.c)*(b.c)*(1 - N.c/(k.c*(1.01-phi))) - ((beta.cc*S.c*I.c)/N.c^kappa + (epsilon*beta.mc*S.c*I.m)/(N.c+epsilon*N.m)^kappa)- d.c*S.c  + sigma.c*R.c
    dI.c<-(beta.cc*S.c*I.c)/N.c^kappa+(epsilon*beta.mc*S.c*I.m)/(N.c+epsilon*N.m)^kappa - I.c*(alpha.c*d.c+d.c+gamma.c)
    dR.c<-I.c*gamma.c-R.c*d.c
    
    dS.m<-N.m*(b.m)*(1-N.m/(k.m*phi)) - ((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.cm*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - d.m*S.m 
    dI.m<-((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.cm*S.m*I.c)/(epsilon*N.c+N.m)^kappa) - I.m*(alpha.m*d.m+d.m+gamma.m)
    dR.m<-I.m*gamma.m-R.m*d.m
    
    return(list(c(dS.c, dI.c, dR.c,dS.m, dI.m, dR.m)))
  })
}

# Plotting Disease dynamics over time
plot.dynamics<- function(x){
  par(mfrow=c(2,2))
  matplot(x[,"time"],x[,2:4],type="l",xlab="time",ylab="number", bty='n',cex=0.8,
          main="Core Hosts",lwd=2, ylim=c(0,100),col=c("black","darkred","forestgreen"))
  legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20,bty='n',cex=0.7)
  matplot(x[,"time"],(x[,3]/(x[,2]+x[,3]+x[,4])),type="l",xlab="time",ylab="number", bty='n',
          main="Prop. Infected Core Hosts",lwd=2, ylim=c(0,1),cex=0.8) 
  matplot(x[,"time"],x[,5:7],type="l",xlab="time",ylab="number", bty='n',
          main="Matrix Hosts",lwd=2, ylim=c(0,100),cex=0.8,col=c("black","darkred","forestgreen"))
  legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20, bty='n',cex=0.7)
  matplot(x[,"time"],(x[,6]/(x[,5]+x[,6]+x[,7])),type="l",xlab="time",ylab="number", bty='n',
          main="Prop. Infected Matrix Hosts",lwd=2, ylim=c(0,1),cex=0.8)
}

