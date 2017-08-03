### cleaning & processing data for plotting

phi = seq(0,1, by=0.01)
n.sims = 1000
epidemic.size<-seq(0,400,length.out=401)
K.c = 400*(1.01-phi)
epsilon = (1+cos(phi*(pi*3/2)-2.3)) 
spillover = K.c * 0.0001 *epsilon


bins = stoch[,2:401]

sum_stoch = as.data.frame(matrix (NA, nrow =length(phi), ncol=4))

for (j in 1:length(phi)){
  vector <- as.vector(NA)
  for (i in 1:length(epidemic.size)){
    v = bins[j,i]*epidemic.size[i]
    vector <- c(vector, v)
    sum_stoch[j,2] = sum(vector, na.rm =T )/n.sims
  }
  sum_stoch[j,3] = (n.sims- bins[j,1])/n.sims
  sum_stoch[j,4] = spillover[j]/sum_stoch[j,2]
}
sum_stoch[,1]= phi
colnames(sum_stoch) = c('phi','mean.epi','p.epi','spill:mean inf')

write.csv(sum_stoch,"stochastic/output/summary_for_3panel.csv", row.names = F)
