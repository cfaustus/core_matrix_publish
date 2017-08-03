library(tgp)

# parameter matrix for core-matrix model where humans are in the matrix habitat
set.seed(101)
s1 <- lhs(1000,rbind(
  #ecological characteristics of patch host
  c(0,5), #rmax
  c(0.05,0.5),  #d.p
  c(100,1000), #k.p
  #disease characteristics of patch host
  c(0.0192,1), #gamma.p
  c(0,10), #alpha.p
  #####################################
  # transmission between patch and matrix
  c(0.01,1), #beta.pm
  c(0.01,1), #beta.mp
  #####################################
  #ecological characteristics of humans
  #c(0,0.1095), #rmax.m
  #c(0.0333,0.0167),  #d.m
  #c(100,1000), #k.m
  #disease characteristics of matrix host
  #c(0.0192,1), #gamma.m
  #c(0,10) #alpha.m
  #c(0.0001,1.0) #forested proportion
)
)

colnames(s1)<-c(
  'rmax.p',
  'd.p',
  'k.p',
  'gamma.p', 
  'alpha.p',
  'beta.pm',
  'beta.mp')
  # 'rmax.m',
  # 'd.m',
  # 'k.m',
  # 'gamma.m',
  # 'alpha.m')

s1 <- as.data.frame(s1)

## FIXED PARAMS
phi=0.5
s1$phi = phi
s1$epsilon= (1+cos(phi*(pi*3/2)-2.3))
s1$kappa=1 # freq dependent
R0_C = 1.5 # R0 in core at f=1.0
R0_M = 0.5 # R0 in matrix at f=1.0

s1$rmax.m <- 0.0547 #c(0,0.1095), #rmax.m
s1$d.m <- 0.0167 #c(0.00333,0.0167),  #d.m
s1$k.m <- 500 #c(100,1000), #k.m
s1$gamma.m <- 0.02 #c(0.0192,1), #gamma.m
s1$alpha.m <- 2 #c(0,10) #alpha.m

## calculate beta.pp and beta.mm based on 100% forest cover
s1$beta.pp <- (R0_C*(s1[['gamma.p']]+s1[['alpha.p']]*s1[['d.p']]+s1[['d.p']]))/s1[['k.p']]
s1$beta.mm <- (R0_M*(s1[['gamma.m']]+s1[['alpha.m']]*s1[['d.m']]))/s1[['k.m']]
#s1$beta.pp <- R0_C*(s1[['gamma.p']]+s1[['alpha.p']]*s1[['d.p']]+s1[['d.p']])
#s1$beta.mm <- R0_M*(s1[['gamma.m']]+s1[['alpha.m']]*s1[['d.m']]+s1[['d.m']])

head(s1)

write.csv(s1,"lhc_nomatrixv_dd.csv", row.names = F)

