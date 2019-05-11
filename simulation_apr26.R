##THIS IS THE SIMULATION##

library(deSolve)
library(parallel)

wd <- getwd()
setwd(wd)
rm(wd)

###############
#CONTORL PANEL#
###############
rm(list=ls())

numCores <- 28
options(mc.cores=numCores)

T <- 99
time <- seq(0,T,1)
d <- 10000

x0 <- 4
y0 <- 1

state <- c(x0 = x0, y0 = y0)

set.seed(7326)
#original parameters
alpha_star <- runif(1,min=0.001,max=5) 
theta_a <- runif(1,min=0.001,max=alpha_star)
beta_star <- runif(1,min=0.001,max=alpha_star)
theta_b <- runif(1,min=0.001,max=beta_star)
delta_star <- runif(1,min=0.001,max=5)
theta_d <- runif(1,min=0.001,max=delta_star)
gamma_star <- runif(1,min=0.001,max=delta_star)
theta_g <- runif(1,min=0.001,max=gamma_star)

parms <- c(alpha_star = alpha_star, theta_a = theta_a, beta_star = beta_star, 
           theta_b = theta_b, delta_star = delta_star, theta_d = theta_d, 
           gamma_star = gamma_star, theta_g = theta_g)

oParm <- parms

#external data
z_t <- round(runif(length(time), min = 1, max = 50), 2)
###############

###########
#FUNCTIONS#
###########
sigimp <- approxfun(time, z_t, rule=2)

#LOTKA VOLTERRA FUNCTION
LotVmod <- function (Time, State, pars) {
  with(as.list(c(State, pars)), {
    alpha <- alpha_star + theta_a*sigimp(Time)
    beta <- beta_star  + theta_b*sigimp(Time)
    gamma <- gamma_star + theta_g*sigimp(Time)
    delta <- delta_star + theta_d*sigimp(Time)
    dx = alpha*State[1]-(beta*State[1]*State[2])
    dy = (-delta*State[2]) + (gamma*State[1]*State[2])
    return(list(c(dx, dy)))
  })
}

#OPTIMIZER FUNCTION
optim_fun <- function(p){
  
  N <- 1
  
  alpha_star <- runif(1,0.001,5)
  beta_star <- runif(1,0.001,alpha_star)
  delta_star <- runif(1,0.001,5)
  gamma_star <- runif(1,0.001,delta_star)
  
  theta_a <- runif(1,0.001,alpha_star)
  theta_b <- runif(1,0.001,beta_star)
  theta_d <- runif(1,0.001,delta_star)
  theta_g <- runif(1,0.001,gamma_star)

  par <- c(alpha_star = alpha_star, theta_a = theta_a, beta_star = beta_star, 
           theta_b = theta_b, delta_star = delta_star, theta_d = theta_d, 
           gamma_star = gamma_star, theta_g = theta_g)
  
  ##CHANGED rtol FROM 1E-15##
  out <- ode(state, time, LotVmod, par, rtol = 1e-15, maxsteps = 500000)
  
  
  out[is.na(out)] <- 0
  
  ssegov <- sum((ts[,5]-out[,3])^2)
  ssebh <- sum((ts[,4] - out[,2])^2) 
  
  sse <- ssegov + ssebh
  
  output <- c(sse = sse, alpha_star = alpha_star, theta_a = theta_a, beta_star = beta_star, 
              theta_b = theta_b, delta_star = delta_star, theta_d = theta_d, 
              gamma_star = gamma_star, theta_g = theta_g)
  
  return(output)
}

print("Functions Processed Successfully")
###########

#################
#DATA GENERATION#
#################
ts <- data.frame(ode(state, time, LotVmod, parms, rtol = 1e-15, maxsteps = 500000))

ts[is.na(data)] <- 0

ts$x <- exp(log(ts$x0)) * exp(rnorm(100, sd=.25))
ts$y <- exp(log(ts$y0)) * exp(rnorm(100, sd=.25))
#################


######################
#OPTIMIZATION ROUTINE#
######################
for(i in seq(1,d,1)){
  fit <- mclapply(1:numCores, optim_fun)
  
  #print(fit)
  
  try(fit1 <- data.frame(fit, stringsAsFactors = FALSE))
  
  try(fit2 <- data.frame(t(fit1)))
  try(rownames(fit2) <- seq(1,nrow(fit2),1))
  try(write.csv(fit2, paste("simulation_output",i,".csv", sep = "")))
}
######################







