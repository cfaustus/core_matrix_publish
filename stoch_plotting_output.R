#plotting stochastic simulations
rm(list = ls())

setwd('/Users/cfaust/Dropbox/Projects/NCEAS/Core-Matrix Interface/code_clean/stochastic')
source('image_scale.R')
library(RColorBrewer)
library(reshape2)
library(plyr)

stoch <- read.csv("output/ddsimulationwithepsilon_phi_400sims_24may.csv", header= T)
cases <- read.csv("output/ddsimulationwithepsilon_phi_400sims_cases_24may.csv")

phi = seq(0,1, by=0.01)
epidemic.size<-seq(0,400,length.out=401)
out.matrix <- as.matrix(stoch[,2:401]) #convert to matrix without f 
par(mar=c(5,5,1,1), mgp=c(3, 1, 0))
fire = heat.colors(100, alpha = 1)
fire = rev(fire)
image(z = log10(out.matrix + 1),
      x = phi, y = epidemic.size, 
      xlab = "proportion habitat converted", ylab = "epidemic size",
      col= fire,
      las = 1)
image.scale(log10(out.matrix + 1), col =fire, horiz = F, las =1, ylim=c(0,3), xaxt='n',xayt='n')
polygon(c(0,0,1), c(400,0, 400), col= 'white', border = NA)


n.sim =10000
phi = seq(0,1, by=0.01)
epidemic.size<-seq(0,400,length.out=401)
K.c = 400*(1.01-phi)
epsilon = (1+cos(phi*(pi*3/2)-2.3)) 
spillover = K.c * 0.0001 *epsilon

p.epi <-(n.sim-stoch[,2])/n.sim # calculating number of outbreaks


total<-t(t(stoch[,2:401])*epidemic.size[1:400])
cumulative.inf <- rowSums(total)
mean.epi <-cumulative.inf/n.sim 

stoch.mat <- as.matrix(stoch[,-1])
dimnames(stoch.mat) = list(rev(phi), epidemic.size[1:400])
stoch.tidy <- melt(stoch.mat)
colnames(stoch.tidy) <- c('phi','epidemicsize','count')
head(stoch.tidy)
stoch.tidy<-na.omit(stoch.tidy)

dat1 <- data.frame(length = rep(stoch.tidy$phi, stoch.tidy$count))
dat2 <- data.frame(length = rep(stoch.tidy$epidemicsize, stoch.tidy$count))
stoch.long <- cbind(dat1,dat2)
colnames(stoch.long) <- c("phi", "epidemic_size")
head(stoch.long)



library(ggplot2)


par(mfrow=c(3,1),
    mgp=c(2.8, 1, 0),
    las=1)
#par(fig=c(0,0.8,0.55,1), new=TRUE)
# probability of an outbreak.
par(mar = c(2,4,0,0))
plot(phi, spillover/mean.epi,
     type='l', pch = 20,
     col = 'navy',
     # bty = 'n',
     ylab= 'spillover to epidemic',
     #ylim= c(0, 0.015),
     xaxt='n',
     xlab = '',
     lwd=3)
#par(fig=c(0,0.8,0.55,1), new=TRUE)
plot(phi, p.epi, 
     type='b', pch = 20,
     col = 'navy',
     # bty = 'n',
     ylab = 'probabilty of matrix epidemic',
     xaxt='n')
#par(fig=c(0.65,1,0,0.8),new=TRUE)
par(mar = c(4,4,0,0))
plot(phi, mean.epi,
     type='b', pch = 20,
     col = 'navy',
     xlab = 'proportion converted habitat', 
     #bty= 'n',
     ylab= 'mean epidemic size',
     las=1)
sum_outbreak = as.data.frame(cbind(phi = phi, p.epi = p.epi, mean.epi = mean.epi, spillover = spillover/mean.epi))
#write.csv(sum_outbreak,"stochastic/output/summary_for_3panel.csv", row.names = F)
#distributions of outbreaks
par(mfrow=c(1,3),
    mgp=c(3.1, 1, 0),
    las=1)
par(mar = c(4,5,0,0))
plot(epidemic.size[1:400], (out.matrix[11,]+1)/n.sim, type= 'p',
     pch=20, bty='n', log='y',
     ylab= 'frequency',
     ylim=c(1/n.sim,1),
     xlab = '',
     col= rgb(red= 184/255,green = 134/255,blue = 11/255, alpha = 0.6))
legend(100,1, 'f = 0.1', bty='n', cex =1.5)
par(mar = c(4,2,0,3))
plot(epidemic.size[1:400], (out.matrix[51,]+1)/n.sim, pch=20,log='y',
     col = rgb(red= 139/255,green = 69/255,blue = 19/255, alpha = 0.6),
     yaxt ='n',bty='n',
     ylim=c(1/n.sim,1),
     xlab = 'outbreak size')
legend(100,1, 'f = 0.5', bty='n', cex =1.5)
par(mar = c(4,1,0,4))
plot(epidemic.size[1:400], (out.matrix[91,]+1)/n.sim, pch=20, log='y',
     col = rgb(red= 85/255,green = 107/255,blue = 47/255, alpha = 0.6),
     ylim=c(1/n.sim,1),
     yaxt='n', bty='n', xlab = '')
legend(100,1, 'f = 0.9', bty='n', cex =1.5)

legend('topright', c('0.9','0.5','0.1'),
       title= 'proportion of core habitat',
       col= c('forestgreen', 'black', 'darkgoldenrod'), pch=20,
       bty = 'n', xlab=''
)



###Violin plots

dat.lmh <- subset(stoch.long, phi %in% c(0.1, 0.5, 0.9))
dat.lmh$phi <- as.factor(dat.lmh$phi)
ggplot(dat.lmh, aes(phi, epidemic_size+1)) + 
  geom_violin(width=1.3, fill = 'darkgoldenrod', colour = NA, alpha=0.9) +
  # scale_fill_manual(
  #   breaks=c(0,1,8),
  #   labels=c('no persistance', '1','invasion')) +
  stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
  stat_summary(fun.y=median, geom="point", size=3, color="black") +
  scale_y_log10() +
  xlab("proportion converted") + 
  ylab("epidemic size") + 
  theme_bw()

dat.range <- subset(stoch.long, stoch.long$phi %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
2
dat.range<- na.omit(dat.range)
ggplot(dat.range, aes(phi, epidemic_size+1)) + 
  geom_violin(width=1.3, fill = 'darkgoldenrod', colour = NA, alpha=0.9) +
  xlab("proportion converted") + 
  ylab("epidemic size + 1") + 
  theme_bw() + 
  stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
  stat_summary(fun.y=median, geom="point", size=3, color="black") +
  scale_y_log10() 

ggplot(dat.range, aes(phi, epidemic_size/400)) + 
  geom_violin(width=1.3, fill = 'darkgoldenrod', colour = NA, alpha=0.9) +
  xlab("proportion converted") + 
  ylab("epidemic size + 1") + 
  theme_bw() + 
  stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
  stat_summary(fun.y=median, geom="point", size=3) 


####
head(stoch.tidy)
stoch.tidy$k <- stoch.tidy$phi*400
ggplot(stoch.tidy, aes(phi, epidemicsize*count/(k*n.sim)))  +
  geom_point(size=0.8,alpha=0.5) +
  stat_summary(fun.y=mean, geom="point", shape=5, size=3,color='darkred') +
  theme_bw()


##

cases$phi <- phi
cases <- cases[-1,]
head(cases)

library(Hmisc)
cases_sum <- transform(cases, mean=rowSums(cases[,1:1000])/n.sim)
p_epi<- 1-rowSums(cases_sum == 0)/n.sim
#<- 1-rowSums(cases_sum < 20)/n.sim
P_epi_mat <- as.data.frame(binconf(x=n.sim-rowSums(cases_sum == 0), n=n.sim))
plot(cases_sum$phi,cases_sum$mean, type = 'l')
plot(phi[-1],P_epi_mat$PointEst, 
     ylim = c(0.9,1), 
     type = 'l',
     lwd = 2, 
     xlab = 'habitat converted',
     ylab = 'probability of outbreak')
lines(phi[-1],P_epi_mat$Lower, lty=3, col= 'darkgrey')
lines(phi[-1],P_epi_mat$Upper, lty=3, col= 'darkgrey')

ddply(cases, c("phi"), summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)))

out_size<- as.data.frame(matrix(NA, nrow = length(phi), ncol = 3))
for i in length(CI(Data$ Fish, 
                   ci=0.95) 
                head(P_epi_mat)
                library(Rmisc)
                for (i in 2:length(phi)){
                  #i = 100
                  data <- as.numeric(cases[i,])
                  a <- mean(data)
                  s <- sd(data)
                  error <- qnorm(0.975)*s/sqrt(n.sim)
                  out_size[i,1] <- a 
                  out_size[i,2] <- a-error
                  out_size[i,3] <- a+error
                }
                colnames(out_size) <- c('mean','lower.ci','upper.ci')
                plot(phi[-1],out_size$mean[-1], 
                     ylim = c(0,100), 
                     type = 'l',
                     lwd = 2, 
                     col = 'navy',
                     xlab = 'habitat converted',
                     ylab = 'number infected')
                lines(phi[-1],out_size$lower.ci[-1], lty=3, col= 'darkgrey')
                lines(phi[-1],out_size$upper.ci[-1], lty=3, col= 'darkgrey')
                