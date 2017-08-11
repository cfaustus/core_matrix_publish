#sensitivity of core-matrix params
#
library(sensitivity)
#out is the dataframe that the sensitivity_func_cluster outputs
out = read.csv("output/out_fd_R0vary_11aug17.csv", header =T)
s1 = read.csv("lhc/lhc_humans_fd_R0vary.csv",header=T) 
var = c("rmax.c", 'd.c', 'k.c', 'gamma.c', 'alpha.c', 'psi', 
         "rmax.m", 'd.m', 'k.m', 'gamma.m', 'alpha.m',
         'beta.c', 'beta.m', 'R0_C', 'R0_M')
s1.var = s1[ , var]
s1.demo = s1.var[,c("rmax.c", 'd.c', 'k.c')]
ratio = pcc(X=s1.var,y=(out['I.c']/sum(out[c('I.c','S.c','R.c')]))/(out['I.m']/sum(out[c('I.m','S.m','R.m')])),
           rank=T, nboot=1e2,conf=0.95)
core = pcc(X=s1.var,y=(out['I.c']/sum(out[c('I.c','S.c','R.c')])),
           rank=T, nboot=1e2,conf=0.95)
matrix = prcc(X=s1.var,y=(out['I.m']/sum(out[c('I.m','S.m','R.m')])),
           rank=T, nboot=1e2,conf=0.95)
ratio.inf = pcc(X=s1.var,y=(out['I.c']/out['I.m']),
           rank=T, nboot=1e4,conf=0.95)
core.inf = pcc(X=s1.var,y=out['I.c'],
           rank=T, nboot=1e4,conf=0.95)
matrix.inf  = pcc(X=s1.var,y=(out['I.m']),
            rank=T, nboot=1e4,conf=0.95) 
gamma.beta  =  pcc(X=s1.demo,y=out['I.m'],
                  rank=F, nboot=1e4,conf=0.95)
# running correlation on number cases vs. prevalence is equivalent 

library(gplots)
barplot2(as.vector(core$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = FALSE, ci.u = core$PRCC[[5]], ci.l = core$PRCC[[4]], 
         col='gray',ci.lwd=3, ci.width = 0, ylim=c(-0.15,0.15), 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='prevalence of matrix infections')

barplot2(as.vector(patch.inf$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = patch$PRCC[[5]], ci.l = patch$PRCC[[4]], 
         col='forestgreen',ci.lwd=3, ci.width = 0, ylim=c(-0.15,0.15), 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='infections in the core')

barplot2(as.vector(matrix.inf$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = matrix.inf$PRCC[[5]], ci.l = matrix.inf$PRCC[[4]], 
         col= 'darkgoldenrod',ci.lwd=3, ci.width = 0, ylim=c(-0.15,0.15), 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='infections in the matrix')


library(gplots)
barplot2(as.vector(patch$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1),
        plot.ci = TRUE, ci.u = patch$PRCC[[5]], ci.l = patch$PRCC[[4]], 
        ci.lwd=2, ci.width = 0, ylim=c(-0.08,0.06), 
        las=2,  cex.names=0.7, ylab=expression(rho),
        main='infections in the core')


####plotting output
head(out,n=4)

par(mfrow=c(2,2))
matplot(out[,"time"],out[,2:4],type="l",xlab="time",ylab="number", bty='n',cex=0.8,
        main="Patch Hosts",lwd=2, ylim=c(0,1000),col=c("black","darkred","forestgreen"))
legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20,bty='n',cex=0.7)
matplot(out[,"time"],(out[,3]/(out[,2]+out[,3]+out[,4])),type="l",xlab="time",ylab="number", bty='n',
        main="Prop. Infected Patch Hosts",lwd=2, ylim=c(0,1),cex=0.8) 
matplot(out[,"time"],out[,5:7],type="l",xlab="time",ylab="number", bty='n',
        main="Matrix Hosts",lwd=2, ylim=c(0,1000),cex=0.8,col=c("black","darkred","forestgreen"))
legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20, bty='n',cex=0.7)
matplot(out[,"time"],(out[,6]/(out[,5]+out[,6]+out[,7])),type="l",xlab="time",ylab="number", bty='n',
        main="Prop. Infected Matrix Hosts",lwd=2, ylim=c(0,1),cex=0.8)




#x = pcc(X=s1,y=(out[,3]/rowSums(out[,2:4]))/(out[,6]/rowSums(out[,5:7])),rank=T)
# #this is how to plot the results
# barplot(as.vector(x$PRCC[[1]]),names.arg = names(s1),cex.names=0.7,las=2,main='Incidences_U/Incidences_D',ylab=expression(rho))
