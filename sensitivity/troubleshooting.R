S.c=(1.01-x[['phi']])*x[['k.c']]*0.7
I.c=1
R.c=0
S.m=x[['phi']]*x[['k.m']]*0.7
I.m=0
R.m=0
N.c = S.c+I.c+R.c
N.m = S.m+I.m+R.m

dS.c <- ifelse((1 - N.c/(x$k.c*(1.01-x$phi))) > 0, #IF UNDER CARRYING CAPACITY THEN..
               (N.c)*(x$rmax.c+x$d.c)*(1 - N.c/(x$k.c*(1.01-x$phi))) 
               - ((x$beta.c*S.c*I.c)/N.c^x$kappa
                  + (x$epsilon*x$psi*x$beta.m*S.c*I.m)/(N.c+x$epsilon*N.m)^x$kappa)
               - x$d.c*S.c 
               + x$sigma.c*R.c, # ELSE... NO BIRTHS
               - x$d.c*S.c 
               + x$sigma.c*R.c) 
dI.c <- (x$beta.c*S.c*I.c)/N.c^x$kappa + (x$epsilon*x$phi*x$beta.m*S.c*I.m)/(N.c+x$epsilon*N.m)^x$kappa  - I.c*(x$alpha.c*x$d.c+x$d.c+x$gamma.c)
dR.c <- I.c*x$gamma.c - R.c*(x$d.c+x$sigma.c)
#Dynamics in the matrix
dS.m <-  ifelse((1 - N.m/(x$k.m*x$phi)) > 0,
                (N.m)*(x$rmax.m+x$d.m)*(1 - N.m/(x$k.m*x$phi)) 
                - ((x$beta.m*S.m*I.m)/N.m^x$kappa
                   + (x$epsilon*x$psi*x$beta.c*S.m*I.c)/(x$epsilon*N.c+N.m)^x$kappa) 
                - x$d.m*S.m 
                + x$sigma.m*R.m, 
                #- ((beta.mm*S.m*I.m)/N.m^kappa
                #    + (epsilon*beta.cm*S.m*I.c)/(epsilon*N.c+N.m)^kappa) 
                - x$d.m*S.m
                + x$sigma.m*R.m)
dI.m <- ((x$beta.m*S.m*I.m)/N.m^x$kappa + (x$epsilon*x$psi*x$beta.c*S.m*I.c)/(x$epsilon*N.c+N.m)^x$kappa) - I.m*(x$alpha.m*x$d.m+x$d.m+x$gamma.m)
dR.m <- I.m*x$gamma.m - R.m*(x$d.m+x$sigma.m)



S.c2 = S.c + dS.c
I.c2 = I.c + dI.c
R.c2 = R.c + dR.c
S.m2 = S.m + dS.m
I.m2 = I.m + dI.m
R.m2 = R.m + dR.m
N.c2 = S.c2+I.c2+R.c2
N.m2 = S.m2+I.m2+R.m2

dS.c2 <- ifelse((1 - N.c2/(x$k.c*(1.01-x$phi))) > 0, #IF UNDER CARRYING CAPACITY THEN..
               ((N.c2)*(x$rmax.c+x$d.c)*(1 - N.c2/(x$k.c*(1.01-x$phi)))- 
                 ((x$beta.c*S.c2*I.c2)/N.c^x$kappa +
                 (x$epsilon*x$psi*x$beta.m*S.c2*I.m2)/(N.c2+x$epsilon*N.m2)^x$kappa) - 
                 x$d.c*S.c2 + x$sigma.c*R.c2), # ELSE... NO BIRTHS
                 (- x$d.c*S.c2 + x$sigma.c*R.c2))
dI.c2 <- (x$beta.c*S.c2*I.c2)/N.c^x$kappa + (x$epsilon*x$phi*x$beta.m*S.c2*I.m2)/(N.c2+x$epsilon*N.m2)^x$kappa  - I.c2*(x$alpha.c*x$d.c+x$d.c+x$gamma.c)
dR.c2 <- I.c2*x$gamma.c - R.c2*(x$d.c+x$sigma.c)
#Dynamics in the matrix
dS.m2 <-  ifelse((1 - N.m2/(x$k.m*x$phi)) > 0,
                (N.m2)*(x$rmax.m+x$d.m)*(1 - N.m2/(x$k.m*x$phi)) 
                - ((x$beta.m*S.m2*I.m2)/N.m2^x$kappa
                   + (x$epsilon*x$psi*x$beta.c*S.m2*I.c2)/(x$epsilon*N.c2+N.m2)^x$kappa) 
                - x$d.m*S.m2
                + x$sigma.m*R.m2, 
                #- ((beta.mm*S.m*I.m2)/N.m2^kappa
                #    + (epsilon*beta.cm*S.m2*I.c)/(epsilon*N.c+N.m2)^kappa) 
                - x$d.m*S.m2
                + x$sigma.m*R.m2)
dI.m2 <- ((x$beta.m*S.m2*I.m2)/N.m2^x$kappa + (x$epsilon*x$psi*x$beta.c*S.m2*I.c2)/(x$epsilon*N.c2+N.m2)^x$kappa) - I.m2*(x$alpha.m*x$d.m+x$d.m+x$gamma.m)
dR.m2 <- I.m2*x$gamma.m - R.m2*(x$d.m+x$sigma.m)



S.c3 = S.c2 + dS.c2
I.c3 = I.c2 + dI.c2
R.c3 = R.c2 + dR.c2
S.m3 = S.m2 + dS.m2
I.m3 = I.m2 + dI.m2
R.m3 = R.m2 + dR.m2
