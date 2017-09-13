# plotting output

plot.dynamics<- function(x,kc,km){
  par(mfrow=c(2,2))
  matplot(x[,"time"],x[,2:4],type="l",xlab="time",ylab="number", bty='n',cex=0.8,
          main="Core Hosts",lwd=2, ylim=c(0,kc),col=c("black","darkred","forestgreen"))
  legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20,bty='n',cex=0.7)
  matplot(x[,"time"],(x[,3]/(x[,2]+x[,3]+x[,4])),type="l",xlab="time",ylab="number", bty='n',
          main="Prop. Infected Core Hosts",lwd=2, ylim=c(0,1),cex=0.8) 
  matplot(x[,"time"],x[,5:7],type="l",xlab="time",ylab="number", bty='n',
          main="Matrix Hosts",lwd=2, ylim=c(0,km),cex=0.8,col=c("black","darkred","forestgreen"))
  legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20, bty='n',cex=0.7)
  matplot(x[,"time"],(x[,6]/(x[,5]+x[,6]+x[,7])),type="l",xlab="time",ylab="number", bty='n',
          main="Prop. Infected Matrix Hosts",lwd=2, ylim=c(0,1),cex=0.8)
}
plot.dynamics(out, x[['k.c']], x[['k.m']])


