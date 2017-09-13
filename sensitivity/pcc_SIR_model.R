#sensitivity of core-matrix params
library(sensitivity)
library(gplots)

out = read.csv("output/out_fd_R0vary_phi0.1_13sept17.csv", header =T)
#out is the dataframe that the sensitivity_func_cluster outputs
s1 = read.csv("sensitivity/lhc/lhc_humans_fd_R0vary_phi0.1.csv",header=T) 
var = c("b.c", 'd.c', 'k.c','R0_C', 'beta.c', 'gamma.c', 'alpha2.c', 
        'psi', 
        "b.m", 'd.m', 'k.m','R0_M', 'beta.m', 'gamma.m', 'alpha2.m')
s1.var = s1[ , var]

# pcc 
bt = 1e2

core = pcc(X=s1.var,y=out['p.c'],
           rank=T, nboot=bt,conf=0.95)
max_core = pcc(X=s1.var,y=out['max.c'],
           rank=T, nboot=bt,conf=0.95)
auc_core = pcc(X=s1.var,y=out['auc.c'],
               rank=T, nboot=bt,conf=0.95)
core.inf = pcc(X=s1.var,y=out['I.c'],
               rank=T, nboot=bt,conf=0.95)
matrix = pcc(X=s1.var,y=out['p.m'],
           rank=T, nboot=bt,conf=0.95)
max_matrix = pcc(X=s1.var,y=out['max.m'],
             rank=T, nboot=bt,conf=0.95)
auc_matrix = pcc(X=s1.var,y=out['auc.m'],
                 rank=T, nboot=bt,conf=0.95)
matrix.inf  = pcc(X=s1.var,y=(out['I.m']),
                  rank=T, nboot=bt,conf=0.95) 

ratio = pcc(X=s1.var,y=out['p.c']/out['p.m'],
            rank=T, nboot=bt,conf=0.95)
ratio.inf = pcc(X=s1.var,y=(out['I.c']/out['I.m']),
           rank=T, nboot=bt,conf=0.95)
# running correlation on number cases vs. prevalence is equivalent 

yrange=c(-0.5,0.8)
tones = c(rep("darkgoldenrod", 7),'grey',rep('forestgreen',7))
par(mfrow=c(2,3))
barplot2(as.vector(core$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = core$PRCC[[5]], ci.l = core$PRCC[[4]], 
         col=tones,ci.lwd=3, ci.width = 0, ylim=yrange, 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='prevalence in the core')
barplot2(as.vector(auc_core$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = auc_core$PRCC[[5]], ci.l = auc_core$PRCC[[4]], 
         col=tones,ci.lwd=3, ci.width = 0, ylim=yrange, 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='auc in the core')
barplot2(as.vector(max_core$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = max_core$PRCC[[5]], ci.l = max_core$PRCC[[4]], 
         col=tones, ci.lwd=3, ci.width = 0, ylim=yrange, 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='max prevalence in the core')
barplot2(as.vector(matrix$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = matrix$PRCC[[5]], ci.l = matrix$PRCC[[4]], 
         col=tones,ci.lwd=3, ci.width = 0, ylim=yrange, 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='prevalence of matrix infections')
barplot2(as.vector(auc_matrix$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = auc_matrix$PRCC[[5]], ci.l = auc_matrix$PRCC[[4]], 
         col=tones,ci.lwd=3, ci.width = 0, ylim=yrange, 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='auc of matrix infections')
barplot2(as.vector(max_matrix$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = max_matrix$PRCC[[5]], ci.l = max_matrix$PRCC[[4]], 
         col=tones, ci.lwd=3, ci.width = 0, ylim=yrange, 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='max prevalence in the matrix')

par(mfrow=c(2,1))
barplot2(as.vector(ratio$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = ratio$PRCC[[5]], ci.l = ratio$PRCC[[4]], 
         col=tones,ci.lwd=3, ci.width = 0, ylim=c(-0.2,0.8), 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='ratio of prevalence')
barplot2(as.vector(ratio.inf$PRCC[[1]]), beside = TRUE, horiz = FALSE, names.arg = names(s1.var),
         plot.ci = TRUE, ci.u = ratio.inf$PRCC[[5]], ci.l = ratio.inf$PRCC[[4]], 
         col=tones,ci.lwd=3, ci.width = 0, ylim=c(-0.2,0.8), 
         las=2,  cex.names=0.7, ylab=expression(rho),
         main='ratio of infections')

