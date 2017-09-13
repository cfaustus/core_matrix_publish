rm(list = ls())

library(snowfall) 
library(parallel)
# detecting cores
ncores <- detectCores()  
#ncores = 4
sfInit(cpus=ncores,parallel=T)
sfSource("sensitivity/sensitivity_func_cluster_phi.R")


out <- sfClusterApplyLB(param.list,super_ode,times=times)
out2 <- as.data.frame(t(sapply(out,function(x) x)))

out3 <- data.frame(matrix(unlist(out2), nrow=1000))
#colnames(out3)<-names(out2)
colnames(out3) = c('S.c','I.c','R.c','S.m', 'I.m','R.m','p.c','p.m','time','auc.c','auc.m','max.c','max.m')
write.csv(out3,"output/out_dd_R0vary_phi0.9_13sept17.csv",row.names = F)

sfStop()
