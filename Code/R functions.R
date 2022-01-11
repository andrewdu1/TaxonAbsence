#######################
#                     #
#     R functions     #
#                     #
#######################

# Author: Andrew Du & Eric Friedlander
# Date: 3-3-21


# function for the standard reflected beta-binomial pdf (eq. S6)
# ARGUMENTS:
# X: number of successes (# Paranthropus specimens)
# n: number of trials (# total mammalIAN specimens)
# lambda: shape parameter
stdReflectBBpdf <- function(X, n, lambda){
  
  require(rmutil)
  
  m <- 1 / (2 - lambda)
  s <- 2 - lambda
  
  # where m*s is the first shape parameter of the beta distribution, and s * (1-m) is the second shape parameter
  
  return(dbetabinom(X, n, m, s))
}

# function for posterior probability that Z_i = 1, i.e., tau (eq. S10)
# ARGUMENTS:
# X: number of successes (# Paranthropus specimens)
# n: number of trials (# total mammalian specimens)
# psi: parameter (unconditional probability that Z_i = 1)
# lambda: standard reflected beta distribution shape parameter
tau <- function(X, n, psi, lambda){
  
  dbetabinom <- stdReflectBBpdf(X, n, lambda)
  
  return((dbetabinom * psi) / ((X == 0) * (1 - psi) + dbetabinom * psi))
}

# function for expected log-likelihood (eq. S11)
# ARGUMENTS:
# same as in tau() function
Q <- function(X, n, psi, lambda){
  
  tau_res <- tau(X, n, psi, lambda)
  
  dbetabinom <- stdReflectBBpdf(X, n, lambda)
  
  return(sum((1 - tau_res) * log(1 - psi) + tau_res * (log(dbetabinom) + log(psi))))
}

# function for estimating lambda in the expected log-likelihood (eq. S11) (to be numerically solved)
# ARGUMENTS:
# X: number of successes (# Paranthropus specimens)
# n: number of trials (# total mammalian specimens)
# psi: parameter (unconditional probability that Z_i = 1)
# lambda: standard reflected beta distribution shape parameter
# tau_res: results from using the tau() function
lambda_Mstep <- function(X, n, psi, lambda, tau_res){
  
  dbetabinom <- stdReflectBBpdf(X, n, lambda)
  
  return(sum(tau_res * log(dbetabinom)))
}

# primary function for the EM algorithm
# ARGUMENTS:
# X: number of successes (# Paranthropus specimens)
# n: number of trials (# total mammalian specimens)
# psi.init: initial guess for psi
# lambda.init: initial guess for lambda
# n.step.max: number of maximum steps for the optimazation process. Only used for our simulations.
EM <- function(X, n, psi.init, lambda.init, n.step.max = NULL){
  
  # step 1: calculate expected log-likelihood given initial guesses of psi and lambda
  Q.res <- Q(X, n, psi.init, lambda.init)
  
  psi.res <- psi.new <- psi.init
  lambda.res <- lambda.new <- lambda.init
  
  delta.param <- 1 # place-holder for changes in parameter estimates after each iteration
  
  while(delta.param > 1e-5){
    
    # step 2: calculate tau using new psi and lambda (E-step)
    tau.res <- tau(X, n, psi.new, lambda.new)
    
    # step 3: estimate new values of psi and lambda given new tau (M-step)
    psi.new <- mean(tau.res)
    
    lambda.opt <- optim(par = lambda.new, fn = lambda_Mstep, gr = NULL, X = X, n = n, psi = psi.new, tau_res = tau.res, method = "L-BFGS-B", control = list(fnscale = -1), lower = -Inf, upper = 0)
    
    lambda.new <- lambda.opt$par[1]
    
    # Calculate revised log-likelihood
    Q.res <- c(Q.res, Q(X, n, psi.new, lambda.new))
    
    # save new parameter estimates
    psi.res <- c(psi.res, psi.new)
    lambda.res <- c(lambda.res, lambda.new)
    
    # Step 4: calculate how much parameters changed
    delta.psi <- psi.res[length(psi.res)] - psi.res[length(psi.res) - 1]
    delta.lambda <- lambda.res[length(lambda.res)] - lambda.res[length(lambda.res) - 1]
    
    # Step 5: stop if difference between old and new parameters is <= 1e-5.
    delta.param <- max(abs(delta.psi), abs(delta.lambda))
    
    # stop while() loop if n.step.max is specified and exceeded
    if(!is.null(n.step.max)) if(length(Q.res) == n.step.max) break
  }
  
  return(list(psi.res = psi.res, p_hat = psi.res[length(psi.res)], lambda.res = lambda.res, lambda_hat = lambda.res[length(lambda.res)], Q.res = Q.res))
}

# function for simulating number of Paranthropus specimens data with known psi and lambda
# ARGUMENTS:
# n_sites: number of sites
# n_mammSpec: vector of number of total mammalian specimens from each site
# psi: psi parameter (probability that a site truly has Paranthropus)
# lambda: shape parameter for beta-binomial pdf
simulateData <- function(n_sites, n_mammSpec, psi, lambda){
  
  require(rmutil) # for simulating random draws from a beta-binomial distribution
  
  Z <- rbinom(n_sites, 1, psi) # Z: whether a site truly has Paranthropus or not
  
  Z_equals_1 <- Z == 1 # index out which sites truly have Paranthropus
  
  X <- Z # save to a new object, X, which is the number of recovered Paranthropus specimens at site i
  
  if(sum(Z_equals_1) > 0) X[Z_equals_1] <- rbetabinom(n = sum(Z_equals_1), size = n_mammSpec[Z_equals_1], m = 1 / (2 - lambda), s = 2 - lambda) # for those sites where Z_i = 1, use a beta-binomial to randomly select a Paranthropus abundance value
  
  res <- data.frame(n_ParanSpec = X, n_mammSpec = n_mammSpec, Z_i = Z) # collate the Paranthropus abundance, total mammal abundance, and Z vectors
  
  return(res)  
}