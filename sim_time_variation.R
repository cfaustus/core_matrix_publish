################################################################################
## script to look at how rate & amount of deforestation affect disease over time
## messy script made by christina.faust@gmail.com 14 april 16
## cleaned...by future christina
################################################################################

rm(list = ls())

source('core_matrix_baseline_func.R')
library(deSolve)
library(plyr)

times <- seq(0, 55, by = 0.01)
times1 <- seq(0, 55, by = 0.01)
times2 <- seq(55.01, 60, by = 0.01)
times3 <- seq(60.01, 65, by = 0.01)
times4 <- seq(65.01, 70, by = 0.01)
times5 <- seq(70.01, 75, by = 0.01)
times6 <- seq(75.01, 80, by = 0.01)
times7 <- seq(80.01, 85, by = 0.01)
times8 <- seq(85.01, 90, by = 0.01)
times9 <- seq(90.01, 100, by = 0.01)

phiseq <- seq(0.1,0.9, by = 0.1)

# CORE
R0.c <- 2.0
gamma.c = 0.01
alpha.c = 0.001
d.c = 0.25
k.c = 200

#MATRIX
R0.m = 0.3
gamma.m = 0.05
alpha.m = 0.01
d.m = 0.05
k.m = 200
phi = phiseq[1]
params.dd <- c(b.c = 0.5, 
               d.c = d.c,
               k.c = k.c,
               beta.cc = R0.c*(gamma.c+alpha.c+d.c), #beta.cc is transmission within core
               beta.cm = 0.2, #beta.cm is transmission from core to matrix
               gamma.c = gamma.c, 
               alpha.c = alpha.c,
               sigma.c = 0.0,
               b.m = 0.1, 
               d.m = d.m,
               k.m = k.m,
               beta.mm = R0.m*(gamma.m+alpha.m+d.m), 
               beta.mc = 0,
               gamma.m = gamma.m, 
               alpha.m = alpha.m,
               sigma.m = 0.0,
               phi = phi,
               epsilon = (1+cos(phi*(pi*3/2)-2.3)),
               kappa = 1)

initial.values <- c(S.c = (1.01-phi)*params.dd[['k.c']]-1,
                    I.c = 1,
                    R.c = 0,
                    S.m = phi*params.dd[['k.m']],
                    I.m = 0,
                    R.m = 0)

out <- as.data.frame(ode(func = patch.matrix.model.phi,
           y = initial.values,
           parms = params.dd,
           times = times1))
#method = 'ode45')
plot.dynamics(out)

##### round 2
max <- length(times)
initial.values2 <- c(S.c = out[max, 'S.c'],
                     I.c = out[max, 'I.c'],
                     R.c = out[max, 'R.c'],
                     S.m = out[max, 'S.m'],
                     I.m = out[max, 'I.m'],
                     R.m = out[max, 'R.m'])
params.dd['phi'] = phiseq[2]
params.dd['epsilon'] = (1+cos(phiseq[2]*(pi*3/2)-2.3))
out2 <- as.data.frame(ode(func = patch.matrix.model.phi,
            y = initial.values2,
            parms = params.dd,
            times = times2)) #, method = 'ode45'

###### round 3
max2 <- length(times2)
initial.values3 <- c(S.c = out2[max2, 'S.c'],
                     I.c = out2[max2, 'I.c'],
                     R.c = out2[max2, 'R.c'],
                     S.m = out2[max2, 'S.m'],
                     I.m = out2[max2, 'I.m'],
                     R.m = out2[max2, 'R.m'])
params.dd['phi'] = phiseq[3]
params.dd['epsilon'] = (1+cos(phiseq[3]*(pi*3/2)-2.3))
out3 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values3,
                          parms = params.dd,
                          times = times3)) #, method = 'ode45'
plot.dynamics(out3)
max3 <- length(times3)
initial.values4 <- c(S.c = out3[max3, 'S.c'],
                     I.c = out3[max3, 'I.c'],
                     R.c = out3[max3, 'R.c'],
                     S.m = out3[max3, 'S.m'],
                     I.m = out3[max3, 'I.m'],
                     R.m = out3[max3, 'R.m'])
params.dd['phi'] = phiseq[4]
params.dd['epsilon'] = (1+cos(phiseq[4]*(pi*3/2)-2.3))
out4 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values4,
                          parms = params.dd,
                          times = times4)) 

max4 <- length(times4)
initial.values5 <- c(S.c = out4[max4, 'S.c'],
                     I.c = out4[max4, 'I.c'],
                     R.c = out4[max4, 'R.c'],
                     S.m = out4[max4, 'S.m'],
                     I.m = out4[max4, 'I.m'],
                     R.m = out4[max4, 'R.m'])
params.dd['phi'] = phiseq[5]
params.dd['epsilon'] = (1+cos(phiseq[5]*(pi*3/2)-2.3))
out5 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values5,
                          parms = params.dd,
                          times = times5)) 

max5 <- length(times5)
initial.values6 <- c(S.c = out5[max5, 'S.c'],
                     I.c = out5[max5, 'I.c'],
                     R.c = out5[max5, 'R.c'],
                     S.m = out5[max5, 'S.m'],
                     I.m = out5[max5, 'I.m'],
                     R.m = out5[max5, 'R.m'])
params.dd['phi'] = phiseq[6]
params.dd['epsilon'] = (1+cos(phiseq[6]*(pi*3/2)-2.3))
out6 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values6,
                          parms = params.dd,
                          times = times6)) 

max6 <- length(times6)
initial.values7 <- c(S.c = out6[max6, 'S.c'],
                     I.c = out6[max6, 'I.c'],
                     R.c = out6[max6, 'R.c'],
                     S.m = out6[max6, 'S.m'],
                     I.m = out6[max6, 'I.m'],
                     R.m = out6[max6, 'R.m'])
params.dd['phi'] = phiseq[7]
params.dd['epsilon'] = (1+cos(phiseq[7]*(pi*3/2)-2.3))
out7 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values7,
                          parms = params.dd,
                          times = times7)) 

max7 <- length(times7)
initial.values8 <- c(S.c = out7[max7, 'S.c'],
                     I.c = out7[max7, 'I.c'],
                     R.c = out7[max7, 'R.c'],
                     S.m = out7[max7, 'S.m'],
                     I.m = out7[max7, 'I.m'],
                     R.m = out7[max7, 'R.m'])
params.dd['phi'] = phiseq[8]
params.dd['epsilon'] = (1+cos(phiseq[8]*(pi*3/2)-2.3))
out8 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values8,
                          parms = params.dd,
                          times = times8)) 

max8 <- length(times8)
initial.values9 <- c(S.c = out8[max8, 'S.c'],
                     I.c = out8[max8, 'I.c'],
                     R.c = out8[max8, 'R.c'],
                     S.m = out8[max8, 'S.m'],
                     I.m = out8[max8, 'I.m'],
                     R.m = out8[max8, 'R.m'])
params.dd['phi'] = phiseq[9]
params.dd['epsilon'] = (1+cos(phiseq[9]*(pi*3/2)-2.3))
out9 <- as.data.frame(ode(func = patch.matrix.model.phi,
                          y = initial.values9,
                          parms = params.dd,
                          times = times9)) 
full<-rbind.fill(out,out2,out3,out4,out5,
                 out6,out7,out8,out9)
plot.dynamics(full)
par(mfrow=c(1,1))
plot(full$time, full$I.c/(full$S.c+full$I.c+full$R.c))
x<-nrow(full)
change <-full[(x/2):x,]
plot(change$time, change$I.c/(change$S.c+change$I.c+change$R.c), 
     ylim=c(0,0.6),
     ylab= 'prevalence',
     type='l',
     col='forestgreen',
     bty='n',lwd=3)
lines(change$time, change$I.m/(change$S.m+change$I.m+change$R.m),
      col='darkgoldenrod', lwd=3)

####POPULATIONS
plot(change$time, (change$S.c+change$I.c+change$R.c), 
     ylim=c(0,150),
     ylab= 'populations',
     type='l',
     col='forestgreen',
     bty='n',lwd=3)
lines(change$time, (change$S.m+change$I.m+change$R.m),
      col='darkgoldenrod', lwd=3)
lines(change$time, 100*change$I.c/(change$S.c+change$I.c+change$R.c), 
      type='l', col='forestgreen',
      bty='n',lwd=3, lty=3)
lines(change$time, 100*change$I.m/(change$S.m+change$I.m+change$R.m),
      col='darkgoldenrod', lwd=3, lty=3)
0
