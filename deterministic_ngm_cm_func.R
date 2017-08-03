# Functions for making F and V matricies for NGM

matrix.F.DD<-function(beta.w, beta.b, epsilon, k, phi){
  z <-matrix(as.numeric(c(beta.w[[1]]*k[[1]]*(1.01-phi), epsilon*beta.b[[1]]*k[[1]]*(1.01-phi),
               epsilon*beta.b[[2]]*k[[2]]*phi, beta.w[[2]]*k[[2]]*phi)),
             nrow = 2, ncol = 2)
  return(z)
}

matrix.F.DD.exdebt<-function(beta.w, beta.b, epsilon, k, phi){
  z <-matrix(c(beta.w[1]*k[1]*sqrt(1.01-phi), epsilon*beta.b[1]*k[1]*sqrt(1.01-phi),
               epsilon*beta.b[2]*k[2]*phi, beta.w[2]*k[2]*phi),
             nrow = 2, ncol = 2)
  return(z)
}
# kx=c(100,100)
# print(test <- matrix.F(params[["beta.w"]],params[["beta.b"]], epsilon[89], kx))
matrix.F.FD<-function(beta.w, beta.b, epsilon){
  z <-matrix(as.numeric(c(beta.w[[1]], epsilon*beta.b[[1]],
               epsilon*beta.b[[2]], beta.w[[2]])),
             nrow = 2, ncol = 2)
  return(z)
}
matrix.V<-function(alpha, gamma, d){
  z <-matrix(c(alpha[1] + gamma[1] + d[1], 0,
               0, alpha[2] + gamma[2] + d[2]), 
             nrow = 2, ncol = 2)
  return(z)
}
# print(test <- matrix.V(params[["alpha"]],params[["gamma"]], params[["d"]]))


plotR0 = function(phi, R0_c, R0_m, R0_combin){
  par(mar = c(5, 5, 3, 3))
  par(mfrow = c(1, 1))
  par(bg = 'white') 
  plot(phi, R0_c, type = "l", lwd = 4, 
       ylim=c(0,2.5),
       bty = 'n', xlim = c(0, 1),
       #main = 'DD with varying between species',
       ylab = expression("R"[0]), xlab = "proportion converted",
       col = 'forestgreen', cex = 1.5)
  lines(phi, R0_m, lwd = 4, col = 'darkgoldenrod')
  abline(h = 1, lty = 2, lwd = 1, col = 'darkred')
  lines(phi, R0_combin[3,], lwd = 2, col = 'gray91')
  lines(phi, R0_combin[5,], lwd = 2, col = 'gray81')
  lines(phi, R0_combin[7,], lwd = 2, col = 'gray71')
  lines(phi, R0_combin[9,], lwd = 2, col = 'gray61')
  lines(phi, R0_combin[11,], lwd = 2, col = 'gray51')
  lines(phi, R0_combin[13,], lwd = 2, col = 'gray41')
  lines(phi, R0_combin[15,], lwd = 2, col = 'gray31')
  lines(phi, R0_combin[17,], lwd = 2, col = 'gray21')
  lines(phi, R0_combin[19,], lwd = 2, col = 'gray11')
  lines(phi, R0_combin[21,], lwd = 2, col = 'gray1')
  
  legend("topleft", 
         title = expression(psi),
         c("1" , "0.5" ,"0.05"),
         lty = c(1, 1, 1), lwd = 3, 
         col = c("gray0", "gray50", "gray90"), 
         bty = 'n' )
}
