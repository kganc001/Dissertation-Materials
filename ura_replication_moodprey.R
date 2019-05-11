##mood/opinion as prey

library(foreign)
library(deSolve)
library(parallel)

###############
#CONTORL PANEL#
###############
rm(list=ls())

start <- Sys.time()
#cmdArgs <- commandArgs(trailingOnly=TRUE)
#numCores <- as.integer(cmdArgs[1])
numCores <- 28
options(mc.cores=numCores)

#number of repetitions
d <- 20000

wd <- getwd()
setwd(wd)
rm(wd)

rep <- read.dta("supreme court mood replication.dta")

#have to transform caselaw because it is negative and can't be negative in LV...
rep$policy_trans <- (rep$policy - min(rep$policy)+1)


#According to thermostatic model (Ura, 110), prey is opinion
ts <- data.frame(rep$mood, rep$policy_trans)

#STATE VARIABLES
x0 <- ts$rep.mood[1] #prey - public opinion/mood
y0 <- ts$rep.policy_trans[1] #predator - policy choice

state <- c(x0 = x0, y0 = y0)

#TIME SERIES
T <- nrow(ts)-1
time <- seq(0,T,1)

#EXTERNAL DATA 
z_t <- data.frame(caselaw = rep$caselaw, unemployment = rep$unemployment, inflation = rep$inflation)
###############


###########
#FUNCTIONS#
###########
sigimp_claw <- approxfun(time, z_t$caselaw, rule=2)
sigimp_unemp <- approxfun(time, z_t$unemployment, rule=2)
sigimp_infl <- approxfun(time, z_t$inflation, rule=2)


#LOTKA VOLTERRA FUNCTION
LotVmod <- function (Time, State, pars) {
  with(as.list(c(State, pars)), {
    alpha <- alpha_star + theta_a_caselaw*sigimp_claw(Time) + theta_a_unemp*sigimp_unemp(Time) + theta_a_infl*sigimp_infl(Time)
    beta <- beta_star  + theta_b_caselaw*sigimp_claw(Time) + theta_b_unemp*sigimp_unemp(Time) + theta_b_infl*sigimp_infl(Time)
    gamma <- gamma_star + theta_g_caselaw*sigimp_claw(Time) + theta_g_unemp*sigimp_unemp(Time) + theta_g_infl*sigimp_infl(Time)
    delta <- delta_star + theta_d_caselaw*sigimp_claw(Time) + theta_d_unemp*sigimp_unemp(Time) + theta_d_infl*sigimp_infl(Time)
    dx = alpha*State[1]-(beta*State[1]*State[2])
    dy = (-delta*State[2]) + (gamma*State[1]*State[2])
    return(list(c(dx, dy)))
  })
}

#OPTIMIZER FUNCTION
optim_fun <- function(p){
  
  N <- 1
  
  alpha_star <- runif(1,0.01,5)
  beta_star <- runif(1,0.01,alpha_star)
  delta_star <- runif(1,0.01,5)
  gamma_star <- runif(1,0.01,delta_star)
  
  theta_a_caselaw <- runif(1,0.01,alpha_star)
  theta_a_unemp <- runif(1,0.01,alpha_star)
  theta_a_infl <- runif(1,0.01,alpha_star)
  
  theta_b_caselaw <- runif(1,0.01,beta_star)
  theta_b_unemp <- runif(1,0.01,beta_star)
  theta_b_infl <- runif(1,0.01,beta_star)
  
  theta_d_caselaw <- runif(1,0.01,delta_star)
  theta_d_unemp <- runif(1,0.01,delta_star)
  theta_d_infl <- runif(1,0.01,delta_star)
  
  theta_g_caselaw <- runif(1,0.01,gamma_star)
  theta_g_unemp <- runif(1,0.01,gamma_star)
  theta_g_infl <- runif(1,0.01,gamma_star)
  
  par <- c(alpha_star = alpha_star,
           theta_a_caselaw = theta_a_caselaw,
           theta_a_unemp = theta_a_unemp,
           theta_a_infl = theta_a_infl,
           beta_star = beta_star,
           theta_b_caselaw = theta_b_caselaw ,
           theta_b_unemp = theta_b_unemp ,
           theta_b_infl = theta_b_infl ,
           delta_star = delta_star ,
           theta_d_caselaw = theta_d_caselaw ,
           theta_d_unemp = theta_d_unemp ,
           theta_d_infl = theta_d_infl ,
           gamma_star = gamma_star ,
           theta_g_caselaw = theta_g_caselaw ,
           theta_g_unemp = theta_g_unemp ,
           theta_g_infl = theta_g_infl)
  
  ##CHANGED rtol FROM 1E-15##
  out <- ode(state, time, LotVmod, par, rtol = 1e-15, maxsteps = 500000)
  
  
  out[is.na(out)] <- 0
  
  ssegov <- sum((ts[,2]-out[,3])^2)
  ssebh <- sum((ts[,1] - out[,2])^2) 
  
  sse <- ssegov + ssebh
  
  output <- c(sse = sse, alpha_star = alpha_star,
              theta_a_caselaw = theta_a_caselaw,
              theta_a_unemp = theta_a_unemp,
              theta_a_infl = theta_a_infl,
              beta_star = beta_star,
              theta_b_caselaw = theta_b_caselaw ,
              theta_b_unemp = theta_b_unemp ,
              theta_b_infl = theta_b_infl ,
              delta_star = delta_star ,
              theta_d_caselaw = theta_d_caselaw ,
              theta_d_unemp = theta_d_unemp ,
              theta_d_infl = theta_d_infl ,
              gamma_star = gamma_star ,
              theta_g_caselaw = theta_g_caselaw ,
              theta_g_unemp = theta_g_unemp ,
              theta_g_infl = theta_g_infl)
  
  return(output)
}

print("Functions Processed Successfully")
###########

######################
#OPTIMIZATION ROUTINE#
######################

for(i in seq(1,d,1)){
  fit <- mclapply(1:numCores, optim_fun)
  
  print(fit)
  
  try(fit1 <- data.frame(fit, stringsAsFactors = FALSE))
  
  try(fit2 <- data.frame(t(fit1)))
  try(rownames(fit2) <- seq(1,nrow(fit2),1))
  try(write.csv(fit2, paste("output_mood",i,".csv", sep = "")))
  
  #  print(paste("Finished Loop ", i))
}

elapsedTime <-  Sys.time() - start
print(elapsedTime)
######################



####################
#FINAL CALCULATIONS#
####################
# minval <- with(data1, data1[sse == min(sse),])
# 
# out <- ode(state, time, LotVmod, minval[1,-1], rtol = 1e-15, maxsteps = 500000)
####################


























