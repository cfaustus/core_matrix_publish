library(tgp)

# parameter matrix for core-matrix model where humans are in the matrix habitat
set.seed(101)
s1 <- lhs(1000,rbind(
  #ecological characteristics of core host
  c(0,5), #rmax
  c(0.05,0.5),  #d.c
  c(100,10000), #k.c
  #disease characteristics of patch host
  c(0.0192,1), #gamma.c
  c(0,10), #alpha.c
  #####################################
  # transmission between patch and matrix
  c(0.01,1), 
  #####################################
  #ecological characteristics of humans
  c(0,0.1095), #rmax.m
  c(0.0333,0.0167),  #d.m
  c(100,1000), #k.m
  #disease characteristics of matrix host
  c(0.0192,1), #gamma.m
  c(0,10) #alpha.m
  #c(0.0001,1.0) #forested proportion
)
)

colnames(s1)<-c(
  'rmax.c',
  'd.c',
  'k.c',
  'gamma.c', 
  'alpha.c',
  'psi',
  'rmax.m',
  'd.m',
  'k.m',
  'gamma.m',
  'alpha.m')

s1 <- as.data.frame(s1)

## FIXED PARAMS
phi=0.5
s1$phi = phi
s1$epsilon= (1+cos(phi*(pi*3/2)-2.3))
s1$kappa=1 # freq dependent
R0_C = 1.5 # R0 in core at phi = 1.0
R0_M = 0.5 # R0 in matrix at phi = 1.0

## calculate beta.cp and beta.mm based on 100% forest cover
#s1$beta.c <- R0_C*(s1[['gamma.c']]+s1[['alpha.c']]*s1[['d.c']]+s1[['d.c']])/(1/s1[['k.c']])
#s1$beta.m <- R0_M*(s1[['gamma.m']]+s1[['alpha.m']]*s1[['d.m']]+s1[['d.m']])/(1/s1[['k.m']]) # is this right?
s1$beta.c <- R0_C*(s1[['gamma.c']]+s1[['alpha.c']]*s1[['d.c']]+s1[['d.c']])
s1$beta.m <- R0_M*(s1[['gamma.m']]+s1[['alpha.m']]*s1[['d.m']]+s1[['d.m']])

head(s1)

write.csv(s1,"sensitivity/lhc/lhc_humans_fd_cR01.5_mR00.5.csv", row.names = F)

