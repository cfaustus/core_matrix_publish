###############################################################################
###     Functions for DD spillover model
###############################################################################

###  First set routine to define one-step function
SIR.onestep <- function (x, params) {
  # state variables defined
  X <- x[2] # susceptible 
  Y <- x[3] # infected
  Z <- x[4] # recovered
  N <- X + Y + Z # total population
  
  # parameters defined
  beta <- params["beta"] # transmission rate defined in params
  mu <- params["mu"]
  gamma <- params["gamma"]
  spill <- params["spill"]
  alpha <- params["alpha"]
  
  ## each individual rate
  rates <- c(birth=mu*N,
             infection=beta*X*Y+spill*X, #infection is a within pop process and spillover
             recovery=gamma*Y,
             sdeath=mu*X,
             ideath=(mu+alpha)*Y,
             rdeath=mu*Z
  )
  
  #### how the numbers in each column change with each event per timestep
  transitions <- list( 
    birth=c(1,0,0),
    infection=c(-1,1,0),
    recovery=c(0,-1,1),
    sdeath=c(-1,0,0),
    ideath=c(0,-1,0),
    rdeath=c(0,0,-1)
  )
  
  ## total event rate
  total.rate <- sum(rates)
  
  ## waiting time (note exponentially distributed random events)
  tau <- rexp(n = 1, rate = total.rate)
  
  ## which event occurs?
  event <- sample.int(n = 6, size = 1, prob = rates/total.rate)
  
  ## number of cases incremented by 1 if event 2 (infection) occurs  
  x + c(tau,transitions[[event]], as.numeric(event == 2))## double square bracket as recursive from list
  
}

###  Simulation for one run

SIR.simul <- function (x, params, maxstep = 10000) {
  output <- array(dim = c(maxstep + 1, 5))
  colnames(output) <- names(x)
  output[1, ] <- x
  k <- 1
  ## loop until either k > maxstep or time >1 (1 year) 
  ## (the loop will COMMENCE unless t>1, so final time >1;
  ## i.e. stops epidemic after one year)
  while ((k <= maxstep) && (x[1] < 1)) {
    k <- k+1
    output[k,] <- x <- SIR.onestep(x,params)
    
  }
  as.data.frame(output[1:k,])
}