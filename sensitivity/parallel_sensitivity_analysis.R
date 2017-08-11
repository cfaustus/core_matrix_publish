library(snowfall) 
library(parallel)
# detecting cores
ncores <- detectCores()  
#ncores = 4
sfInit(cpus=ncores,parallel=T)
sfSource("sensitivity_func_cluster_phi.R")
out <- sfClusterApplyLB(param.list,super_ode,times=times)
out2 <- as.data.frame(t(sapply(out,function(x) x)))
out3 <- data.frame(matrix(unlist(out2), nrow=1000, byrow=T))
colnames(out3)<-names(out2)
write.csv(out3,"output/out_fd_R01.5core_R00.5matrix_11aug17.csv",row.names = F)

sfStop()
