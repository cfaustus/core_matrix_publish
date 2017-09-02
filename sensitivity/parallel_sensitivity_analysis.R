library(snowfall) 
library(parallel)
# detecting cores
ncores <- detectCores()  
#ncores = 4
sfInit(cpus=ncores,parallel=T)
sfSource("sensitivity/sensitivity_func_cluster_phi.R")
s1<- read.csv("sensitivity/lhc/lhc_humans_fd_R0vary_phi0.1.csv",header=T) # dataframe of random parameters
param.list <- apply(s1,1,as.list) #converting into a list
times <- seq(0, 100, by = 0.1) #time steps are constant for each run

out <- sfClusterApplyLB(param.list,super_ode,times=times)
out2 <- as.data.frame(t(sapply(out,function(x) x)))
#out3 <- data.frame(matrix(unlist(out2), nrow=1000, byrow=T))
#colnames(out3)<-names(out2)
colnames(out2) = 
write.csv(out2,"output/out_fd_R0vary_phi0.1_1sept17.csv",row.names = F)

sfStop()
