## Figures for publication
## finalised may 2017
## christina.faust@gmail.com
rm(list=ls()) # clear workspace
setwd("~/Dropbox/Projects/NCEAS/Core-Matrix Interface/code_clean")

## libraries


############################################################
## figure 3: R0
## 2 panels
phi = seq(0.01, 1, by = 0.001) #vector for different converted proportions
p1 =read.csv("deterministic/output/scenario1_DD.csv")
p2 =read.csv("deterministic/output/scenario2_DD.csv")
p3 =read.csv("deterministic/output/scenario3_DD.csv")
p4 =read.csv("deterministic/output/scenario1_FD.csv")
p5 =read.csv("deterministic/output/scenario2_FD.csv")
p6 =read.csv("deterministic/output/scenario3_FD.csv")
pt_size = 1.2
axis_size = 1.2
plotR0_inner = function(phi, datf){
  #par(mar = c(5, 5, 3, 3))
  par(bg = 'white') 
  plot(phi, datf[22,], type = "l", lwd = 3, 
       ylim=c(0,3.0),
       xaxt ='n', las =1,
       bty = 'o', xlim = c(0, 1),
       #main = 'DD with varying between species',
       ylab = expression("R"[0]), xlab = " ",
       col = 'forestgreen', 
       cex.axis = axis_size, cex = pt_size, cex.lab =axis_size)
  axis(1, seq(0,1,by=0.1), labels =F)
  lines(phi, datf[23,], lwd = 3, col = 'darkgoldenrod')
  lines(phi, datf[3,], lwd = 2, col = 'gray91')
  lines(phi, datf[5,], lwd = 2, col = 'gray81')
  lines(phi, datf[7,], lwd = 2, col = 'gray71')
  lines(phi, datf[9,], lwd = 2, col = 'gray61')
  lines(phi, datf[11,], lwd = 2, col = 'gray51')
  lines(phi, datf[13,], lwd = 2, col = 'gray41')
  lines(phi, datf[15,], lwd = 2, col = 'gray31')
  lines(phi, datf[17,], lwd = 2, col = 'gray21')
  lines(phi, datf[19,], lwd = 2, col = 'gray11')
  lines(phi, datf[21,], lwd = 2, col = 'gray1')
  abline(h = 1, lty = 2, lwd = 1.5, col = 'darkred')
}
plotR0_bottom = function(phi, datf){
  #par(mar = c(5, 5, 3, 3))
  par(bg = 'white') 
  plot(phi, datf[22,], type = "l", lwd = 3, 
       ylim=c(0,3.0),
       xaxt ='n', las = 1,
       bty = 'o', xlim = c(0, 1),
       #main = 'DD with varying between species',
       xlab = expression('proportion habitat converted ' (phi) ),
       ylab = expression("R"[0]), 
       col = 'forestgreen', 
       cex.axis = axis_size, cex = pt_size, cex.lab =axis_size)
  axis(1, seq(0,1,by=0.1), cex.axis = axis_size)
  mtext(expression('proportion habitat converted ' (phi) ),
        side=1,line=2.5, cex =axis_size*0.8)
  lines(phi, datf[23,], lwd = 3, col = 'darkgoldenrod')
  lines(phi, datf[3,], lwd = 2, col = 'gray91')
  lines(phi, datf[5,], lwd = 2, col = 'gray81')
  lines(phi, datf[7,], lwd = 2, col = 'gray71')
  lines(phi, datf[9,], lwd = 2, col = 'gray61')
  lines(phi, datf[11,], lwd = 2, col = 'gray51')
  lines(phi, datf[13,], lwd = 2, col = 'gray41')
  lines(phi, datf[15,], lwd = 2, col = 'gray31')
  lines(phi, datf[17,], lwd = 2, col = 'gray21')
  lines(phi, datf[19,], lwd = 2, col = 'gray11')
  lines(phi, datf[21,], lwd = 2, col = 'gray1')
  abline(h = 1, lty = 2, lwd = 1.5, col = 'darkred')
}

f3 <- rbind(c(7,8),
          c(1, 4),c(1, 4),c(1, 4),c(1, 4),
          c(2, 5),c(2, 5),c(2, 5),c(2, 5),
          c(3, 6),c(3, 6),c(3, 6),c(3, 6))
layout(f3)
layout.show(8)

par(omi=c(0.5,0.3,0,0), plt=c(0.15,0.9,0.0,0.9))
#par(mar = c(0.5,6,0.5,1))
#par(mar = c(0.5,0,0.5,0))
plotR0_inner(phi,p1[,-1])
legend("topleft", 
       title = expression('intraspecies efficiency, ' ~ psi),
       c("1" , "0.5" ,"0.05"),
       lty = c(1, 1, 1), lwd = 3,
       col = c("gray0", "gray50", "gray90"),
       bty = 'n', cex = 1.2 )
text(.98,2.9, labels = "A", cex = 2)
plotR0_inner(phi,p2[,-1])
text(.98,2.9, labels = "B", cex = 2)
plotR0_bottom(phi,p3[,-1])
text(.98,2.9, labels = "C", cex = 2)
plotR0_inner(phi,p4[,-1])
legend("topleft", 
       title = expression('single host R'[0]),
       c("core species", "matrix species"),
       lty = c(1, 1), lwd = 3,
       col = c("forestgreen", "darkgoldenrod"),
       bty = 'n', cex = 1.2 )
text(.98,2.9, labels = "D", cex = 2)
plotR0_inner(phi,p5[,-1])
text(.98,2.9, labels = "E", cex = 2)
plotR0_bottom(phi,p6[,-1])
text(.98,2.9, labels = "F", cex = 2)

#par(mar = c(.5,.5,.5,.5), plt=c(1,1,1,.1),omi=c(0,0,0,0))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
text(1, -0.1, labels = "density dependent transmission", cex = 1.5)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
text(1, -0.1, labels = "frequency dependent transmission", cex = 1.5)


############################################################
## figure 4: stochastic simulations

f4 <- rbind(c(1, 1, 1, 4, 4, 4, 4, 4, 5),
           c(2, 2, 2, 4, 4, 4, 4, 4, 5),
           c(3, 3, 3, 4, 4, 4, 4, 4, 5))
layout(f4)
layout.show(5)
phi = seq(0,1.0, by =0.1)
par(mar = c(5,6,1,1))

axis_size = 1.5
pt_size = 1.2
s_3 = read.csv("stochastic/output/summary_for_3panel.csv", header = T)
head(s_3)
plot(s_3$phi,s_3$spill.mean.inf, xaxt ='n', xlab = '',
     bty = 'o', ylab ='spillover : matrix infections',
     pch = 19, type = 'b', col = 'navy',las = 1,
     cex.axis = axis_size, cex = pt_size, cex.lab =axis_size) 
text(0.97,0.8, labels = "A", cex = 2)
axis(1, seq(0,1,by=0.1), labels =F)
#axis(2,cex.axis=1.2)
#mtext("Awesome Y variable", side=2, line=2.2, cex=2, las =1)
plot(s_3$phi,s_3$p.epi, xaxt ='n', xlab = '',
     pch = 19, type = 'b', col = 'navy',las = 1,
     bty = 'o', ylab ='probability of matrix infections',
     cex.axis = axis_size, cex = pt_size, cex.lab =axis_size)
text(0.97,.95, labels = "B", cex = 2)
axis(1, seq(0,1,by=0.1), labels =F)
plot(s_3$phi,s_3$mean.epi, 
     bty = 'o', ylab ='mean epidemic size',
     pch = 19, type = 'b', col = 'navy',las = 1,
     xlab = expression('proportion habitat converted ' (phi) ),
     cex.axis = axis_size, cex = pt_size, cex.lab =axis_size)
text(0.97,90, labels = "C", cex = 2)
axis(1, seq(0,1,by=0.1), labels =F)
source('stochastic/image_scale.R')
stoch <- read.csv("stochastic/output/ddsimulationwithepsilon_phi_400sims_24may.csv", header= T)
phi = seq(0,1, by=0.01)
epidemic.size<-seq(0,400,length.out=401)
out.matrix <- as.matrix(stoch[,2:401]) #convert to matrix without f 
par(mar=c(5,5,1,1), mgp=c(3, 1, 0))
fire = heat.colors(100, alpha = 1)
fire = rev(fire)
image(z = log10(out.matrix + 1),
      x = phi, y = epidemic.size, 
      ylab = "epidemic size",
      col= fire,
      las = 1,
      xlab = expression('proportion habitat converted ' (phi) ),
      cex.axis = axis_size, cex = pt_size, cex.lab =axis_size)
text(.95,390, labels = "D", cex = 2)
axis(1, seq(0,1,by=0.1), labels =F)
par(mar = c(12,4,12,4))
image.scale(log10(out.matrix + 1), col =fire, horiz = F, las =1, ylim=c(0,3), xaxt='n',xayt='n')
mtext("number of simulations", line=1.5, side=2)
mtext("0", line=0.5, side=1)
mtext("1000", line=0.5, side=3)

############################################################
## figure 5: violin plot
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

############################################################
## figure 6: time variations
