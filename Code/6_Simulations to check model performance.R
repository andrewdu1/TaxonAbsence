##################################################
#                                                #
#     Simulations to check model performance     #
#                                                #
##################################################

# Author: Andrew Du


# Read in dataset
d <- read.csv("Datasets/NISP data.csv", header = TRUE)

# Source functions
source("Code/1_R functions.R")

# define objects
x <- d$Paran_nisp
n <- x + d$NonParanMamm_nisp


# set simulated parameter values
psi.sim <- seq(0.1, 0.9, 0.2)
lambda.sim <- c(1, 3, 5, 10, 25, 50, 100, 200)
alpha.sim <- seq(0.1, 1.1, 0.2)
beta.sim <- c(10, 25, 50, 100, 250, 500, 1000)

# number of iterations for our simulations
n.iter <- 1000


# need to do parallel computing for the simulations
# load packages
library(doParallel)
library(parallel)

# get number of cores on computer
numCores <- detectCores()

# make the cluster using all cores
cl <- makeCluster(numCores)

# register the parallel backend
registerDoParallel(cl)

# Two-parameter zero-inflated beta-binomial model
set.seed(100)

BB2.sim.res <- foreach(psi.iter = psi.sim) %:% # iterate through p
  
  foreach(lambda.iter = lambda.sim) %:% # iterate through lambda
  
  foreach(i = seq_len(n.iter), .combine = "rbind") %dopar% { # run through each iteration
    
    x.sim <- simulateData(
      n_sites = nrow(d), 
      n_mammSpec = n, 
      psi = psi.iter, 
      beta_param = lambda.iter) # simulate Paranthropus abundance data for each site, using the observed number of mammalian specimens at each site (n_i) & known values of psi and lambda
    
    while(sum(x.sim$n_ParanSpec) == 0){ 
      
      x.sim <- simulateData(
        n_sites = nrow(d), 
        n_mammSpec = n, 
        psi = psi.iter, 
        beta_param = lambda.iter)
    } # in case some simulations have zero Paranthropus NISP across all sites 

    ZI.BB2.fit <- 
      try(
        ZI.BB2.mle(
          x = x.sim$n_ParanSpec, 
          n = n, 
          psi.start = 0.5,
          lambda.start = 100), 
        silent = TRUE) # try() identifies whether model fitting results in an error (i.e., non-finite log-likelihood)
    
    if(class(ZI.BB2.fit) == "try-error"){
      
      return(c(NA, NA)) # if error, return NAs
      
    } else {
      
      return(ZI.BB2.fit$par) # return vector of estimated parameters
    }
  }

# turn off cluster
stopCluster(cl)

# wrangle data in usable format
BB2.sim.res1 <- array(
  dim = c(
    length(psi.sim), 
    length(lambda.sim), 
    n.iter, 
    2), # psi & lambda estimates
  dimnames = list(
    psi.sim,
    lambda.sim,
    seq_len(n.iter),
    c("psi_hat", "lambda_hat")
  ))

for(i in seq_along(psi.sim)){
  for(j in seq_along(lambda.sim)){
    
    BB2.sim.res1[i, j, , ] <- BB2.sim.res[[i]][[j]]
  }
}

# save the simulation results.
#saveRDS(BB2.sim.res1, file = "Results/BB2 sim results.rds")


# register the parallel backend
registerDoParallel(cl)

# Three-parameter zero-inflated beta-binomial model
set.seed(100)

BB3.sim.res <- foreach(psi.iter = psi.sim) %:% # iterate through psi
  
  foreach(alpha.iter = alpha.sim) %:% # iterate through alpha
  
  foreach(beta.iter = beta.sim) %:% # iterate through beta
  
  foreach(i = seq_len(n.iter), .combine = "rbind") %dopar% { # run through each iteration
    
    x.sim <- simulateData(
      n_sites = nrow(d), 
      n_mammSpec = n, 
      psi = psi.iter, 
      beta_param = c(alpha.iter, beta.iter)) # simulate Paranthropus abundance data for each of our 56 sites, using the observed number of mammalian specimens at each site (n_i) and the known values of psi and lambda
    
    while(sum(x.sim$n_ParanSpec) == 0){ 
      
      x.sim <- simulateData(
        n_sites = nrow(d), 
        n_mammSpec = n, 
        psi = psi.iter, 
        beta_param = c(alpha.iter, beta.iter))
    } # in case some simulations have zero Paranthropus NISP across all sites 
    
    ZI.BB3.fit <- 
      try(
        ZI.BB3.mle(
          x = x.sim$n_ParanSpec,
          n = n,
          psi.start = 0.5,
          alpha.start = 1,
          beta.start = 100
        ), 
        silent = TRUE) # try() identifies whether model fitting results in an error (i.e., non-finite log-likelihood)
    
    if(class(ZI.BB3.fit) == "try-error"){
      
      return(c(NA, NA, NA)) # if error, return NAs
      
    } else {
      
      return(ZI.BB3.fit$par) # return vector of estimated parameters
    }
  }

# turn off cluster
stopCluster(cl)

# wrangle data in usable format
BB3.sim.res1 <- array(
  dim = c(
    length(psi.sim), 
    length(alpha.sim), 
    length(beta.sim),
    n.iter, 
    3), # psi, alpha, & beta estimates
  dimnames = list(
    psi.sim,
    alpha.sim,
    beta.sim,
    seq_len(n.iter),
    c("psi_hat", "alpha_hat", "beta_hat")
  ))

for(i in seq_along(psi.sim)){
  for(j in seq_along(alpha.sim)){
    for(k in seq_along(beta.sim)){
      
      BB3.sim.res1[i, j, k, , ] <- BB3.sim.res[[i]][[j]][[k]]
    }
  }
}

# save the simulation results.
#saveRDS(BB3.sim.res1, file = "Results/BB3 sim results.rds")