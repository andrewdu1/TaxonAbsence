#######################
#                     #
#     R functions     #
#                     #
#######################

# Author: Andrew Du & Eric Friedlander


# function for the beta-binomial pmf (eq. SX)
  # ARGUMENTS:
    # x: number of successes (# Paranthropus specimens)
    # n: number of trials (# total mammalian specimens)
    # lambda: shape parameter (must be >=0)
BBpmf <- function(x, n, lambda){
  
  require(rmutil)
  
  m <- 1 / (2 + lambda)
  s <- 2 + lambda
  
  # where m*s is the first shape parameter of the beta distribution,
  # and s * (1-m) is the second shape parameter
  
  return(dbetabinom(x, n, m, s))
}


# function for calculating log-likelihood of 
# zero-inflated beta-binomial model (eq. SX)
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # param: two-element vector, where the first is psi and the second is lambda
ZI.BB.logL <- function(x, n, param){
  
  psi <- param[1]
  lambda <- param[2]
  
  f_B <- BBpmf(x, n, lambda)
  
  return(sum(log((x == 0) * (1 - psi) + f_B * psi)))
}


# function for computing MLE of psi & lambda
  # ARGUMENTS:
    # psi.start: initial value for psi 
    # lambda.start: initial value for lambda
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # hessian: if TRUE, returns Hessian matrix
ZI.BB.mle <- function(psi.start, lambda.start, x, n, hessian = FALSE){
  
  mle.res <- optim(par = c(psi.start, lambda.start), 
                   fn = ZI.BB.logL,
                   method = "L-BFGS-B",
                   x = x,
                   n = n,
                   lower = c(1e-6, 0),
                   upper = c(1, Inf),
                   control = list(
                     fnscale = -1 # maximize function
                   ),
                   hessian = hessian)
  
  param.hat <- mle.res$par
  names(param.hat) <- c("psi", "lambda")
  
  if(hessian){
    return(list(
      param = param.hat, 
      hessian = mle.res$hessian))
  } else {
    return(param.hat)
  }
}


# function for calculating posterior probability (eq. SX)
# P(Z_i=0 | x_i=0, n_i, psi, lambda)
  # ARGUMENTS:
    # n: number of trials (# total mammalian specimens)
    # psi: psi value
    # lambda: lambda value
post_prob <- function(n, psi, lambda){
  
  return((1 - psi) / (1 - psi + psi * ((1 + lambda) / (n + 1 + lambda))))
}


# Function for bootstrapping data for calculating CIs
# of parameter estimates and posterior probabilities.
# Returns results as a list with n.iter elements
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # n.iter: number of bootstrap iterations
    # silent: if FALSE, prints iteration number
boot_data <- function(x, n, n.iter = 1000, silent = FALSE){
  
  boot.res <- vector(mode = "list", length = n.iter) # empty list to save results to
  
  for(i in seq_along(boot.res)){
    
    index <- sample.int(
      length(x), 
      size = length(x), 
      replace = TRUE) # randomly choose index by which to resample sites
    
    x_boot <- x[index] # get resampled x
    n_boot <- n[index] # get resampled n
    
    x_boot_spec <- numeric(length(x)) # empty vector to save bootstrapped specimens to
    
    for(j in seq_along(x_boot_spec)) x_boot_spec[j] <- 
      rbinom(
        n = 1, 
        size = n_boot[j], 
        prob = x_boot[j] / n_boot[j]) # for each bootstrapped x and n, resample x (number of specimens)
    
    boot.res[[i]] <- data.frame( # save bootstrapped results
      x = x_boot_spec,
      n = n_boot
    )
    
    if(!silent) print(i)
  }
  return(boot.res)
}


# function for simulating number of Paranthropus 
# specimens with known psi and lambda
  # ARGUMENTS:
    # n_sites: number of sites
    # n_mammSpec: vector of number of total mammalian specimens from each site
    # psi: psi value
    # lambda: lambda value
simulateData <- function(n_sites, n_mammSpec, psi, lambda){
  
  require(rmutil) # for simulating random draws from a beta-binomial distribution
  
  Z <- rbinom(n_sites, 1, psi) # Z: whether a site truly has Paranthropus or not
  
  Z_equals_1 <- Z == 1 # index out which sites truly have Paranthropus
  
  X <- Z # save to a new object, X, which is the number of recovered Paranthropus specimens at site i
  
  if(sum(Z_equals_1) > 0){
    X[Z_equals_1] <- rbetabinom(
      n = sum(Z_equals_1), 
      size = n_mammSpec[Z_equals_1], 
      m = 1 / (2 + lambda), 
      s = 2 + lambda) # for those sites where Z_i = 1, use a beta-binomial to randomly select a Paranthropus abundance value
  }
  
  res <- data.frame(
    n_ParanSpec = X, 
    n_mammSpec = n_mammSpec, 
    Z_i = Z) # collate the Paranthropus abundance, total mammal abundance, and Z vectors
  
  return(res)  
}