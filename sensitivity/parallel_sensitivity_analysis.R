library(snowfall) 
# detecting cores
#numCores <- detectCores()  
numCores=8
sfInit(cpus=numCores,parallel=T)
sfSource("sensitivity_func_cluster_phi.R")
out <- sfClusterApplyLB(param.list,super_ode,times=times)
out2 <- as.data.frame(t(sapply(out,function(x) x)))
out3 <- data.frame(matrix(unlist(out2), nrow=1000, byrow=T))
colnames(out3)<-names(out2)
write.csv(out3,"out_fd_nochangeinmatrix_13april16.csv",row.names = F)

sfStop()