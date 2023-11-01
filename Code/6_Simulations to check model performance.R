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

# Define three-parameter beta-binomial model
  # Argument param is a vector with three elements
  # First element is psi
  # Second & third elements are beta-binomial parameters (alpha & beta)
ZI.BB2.logL <- function(x, n, param){

  require(rmutil)
  
  psi <- param[1]
  m <- param[2] / (param[2] + param[3])
  s <- param[2] + param[3]
  
  return(sum(log((x == 0) * (1 - psi) + dbetabinom(x, n, m, s) * psi)))
}

# define objects
x <- d$Paran_nisp
n <- x + d$NonParanMamm_nisp

# create function for simulating data
# function for simulating number of Paranthropus 
# specimens with known psi and lambda
  # ARGUMENTS:
    # n_sites: number of sites
    # n_mammSpec: vector of number of total mammalian specimens from each site
    # psi: psi value
    # A: alpha value (for beta-binomial)
    # B: beta value (for beta-binomial)
simulateData <- function(n_sites, n_mammSpec, psi, A, B){
  
  require(rmutil) # for simulating random draws from a beta-binomial distribution
  
  Z <- rbinom(n_sites, 1, psi) # Z: whether a site truly has Paranthropus or not
  
  Z_equals_1 <- Z == 1 # index out which sites truly have Paranthropus
  
  X <- Z # save to a new object, X, which is the number of recovered Paranthropus specimens at site i
  
  if(sum(Z_equals_1) > 0){
    X[Z_equals_1] <- rbetabinom(
      n = sum(Z_equals_1), 
      size = n_mammSpec[Z_equals_1], 
      m = A / (A + B), 
      s = A + B) # for those sites where Z_i = 1, use a beta-binomial to randomly select a Paranthropus abundance value
  }
  
  res <- data.frame(
    n_ParanSpec = X, 
    n_mammSpec = n_mammSpec, 
    Z_i = Z) # collate the Paranthropus abundance, total mammal abundance, and Z vectors
  
  return(res)  
}

# set simulated parameter values
psi.sim <- seq(0.1, 0.9, 0.2)
alpha.sim <- seq(0.1, 0.9, 0.2)
beta.sim <- seq(10, 1000, length.out = 10)

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

# set seed to make simulations replicable
set.seed(100)

sim.res <- foreach(psi.iter = psi.sim) %:% # iterate through psi
  
  foreach(alpha.iter = alpha.sim) %:% # iterate through alpha
  
  foreach(beta.iter = beta.sim) %:% # iterate through beta
  
  foreach(i = seq_len(n.iter), .combine = "rbind") %dopar% { # run through each iteration
    
    X.sim <- simulateData(
      n_sites = nrow(d), 
      n_mammSpec = n, 
      psi = psi.iter, 
      A = alpha.iter,
      B = beta.iter) # simulate Paranthropus abundance data for each of our 56 sites, using the observed number of mammalian specimens at each site (n_i) and the known values of psi and lambda
    
    ZI.BB2.fit <- optim(par = c(0.5, 0.5, 100), 
                        fn = ZI.BB2.logL,
                        method = "L-BFGS-B",
                        x = X.sim$n_ParanSpec,
                        n = n,
                        lower = c(1e-6, 1e-6, 1e-6),
                        upper = c(1, 1e10, 1e10),
                        control = list(
                          fnscale = -1
                        )) # use MLE to estimate parameters
    
    return(ZI.BB2.fit$par) # return vector of estimated parameters
  }

# turn off cluster
stopCluster(cl)

# save the simulation results.
save(sim.res, file = "Results/Simulation results.rds")