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
ZI.BB2.logL <- function(x, n, param){
  
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
ZI.BB2.mle <- function(psi.start, lambda.start, x, n, hessian = FALSE){
  
  mle.res <- optim(par = c(psi.start, lambda.start), 
                   fn = ZI.BB2.logL,
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
    # beta_param: a vector for parameters of beta distribution. 
      # If length(beta_param) = 1, parameter is lambda.
      # If length(beta_param) = 2, parameters are alpha & beta of the beta distribution
simulateData <- function(n_sites, n_mammSpec, psi, beta_param){
  
  require(rmutil) # for simulating random draws from a beta-binomial distribution
  
  Z <- rbinom(n_sites, 1, psi) # Z: whether a site truly has Paranthropus or not
  
  Z_equals_1 <- Z == 1 # index out which sites truly have Paranthropus
  
  X <- Z # save to a new object, X, which is the number of recovered Paranthropus specimens at site i
  
  # define parameterization of beta-binomial
  if(length(beta_param) == 1){ # one-parameter model
    s <- 2 + beta_param
    m <- 1 / s
  }
  
  if(length(beta_param) == 2){ # two-parameter model
    s <- sum(beta_param)
    m <- beta_param[1] / s
  }
  
  if(sum(Z_equals_1) > 0){
    X[Z_equals_1] <- rbetabinom(
      n = sum(Z_equals_1), 
      size = n_mammSpec[Z_equals_1], 
      m = m, 
      s = s) # for those sites where Z_i = 1, use a beta-binomial to randomly select a Paranthropus abundance value
  }
  
  res <- data.frame(
    n_ParanSpec = X, 
    n_mammSpec = n_mammSpec, 
    Z_i = Z) # collate the Paranthropus abundance, total mammal abundance, and Z vectors
  
  return(res)  
}



########## functions for other models and model selection ########## 
# Binomial model log-likelihood
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # p: probability of success parameter
binom.logL <- function(x, n, p){
  
  return(sum(dbinom(x, n, p, log = TRUE)))
}

# fit binomial model via MLE
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # p.start: initial guess of p to start from
binom.mle <- function(x, n, p.start){
  return(optim(par = p.start, 
               fn = binom.logL,
               method = "L-BFGS-B",
               x = x,
               n = n,
               lower = 1e-6,
               upper = 1 - 1e-6,
               control = list(
                fnscale = -1
               )))
}


# one-parameter beta-binomial model log-likelihood
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # lambda: single parameter of beta-binomial
BB.logL <- function(x, n, lambda){
  
  require(rmutil)
  
  m <- 1 / (2 + lambda)
  s <- 2 + lambda
  
  return(sum(dbetabinom(x, n, m, s, log = TRUE)))
}

# fit one-parameter beta-binomial model via MLE
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # lambda.start: initial guess of lambda to start from
BB.mle <- function(x, n, lambda.start){
  return(optim(par = lambda.start, 
               fn = BB.logL,
               method = "L-BFGS-B",
               x = x,
               n = n,
               lower = 0,
               upper = 1e6,
               control = list(
                 fnscale = -1
               )))
}


# zero-inflated binomial model log-likelihood
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens) 
    # param: two-element vector c(psi, p)
ZI.binom.logL <- function(x, n, param){
  
  psi <- param[1]
  p <- param[2]
  
  return(sum(log((x == 0) * (1 - psi) + dbinom(x, n, p) * psi)))
}

# fit zero-inflated binomial model via MLE
# NB: added in the option of changing optim() parameter search bounds 
# b/c some sites have very low observed relative abundances. Otherwise,
# optim() fails to converge
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # psi.start: initial guess for psi
    # p.start: initial guess for p
    # psi.bounds: bounds for searching for psi values
    # p.bounds: bounds for searching for p values
ZI.binom.mle <- function(x, n, psi.start, p.start, psi.bounds, p.bounds){
  
  return(optim(par = c(psi.start, p.start), 
               fn = ZI.binom.logL,
               method = "L-BFGS-B",
               x = x,
               n = n,
               lower = c(psi.bounds[1], p.bounds[1]),
               upper = c(psi.bounds[2], p.bounds[2]), 
               control = list(
                 fnscale = -1,
                 ndeps = c(1e-4, 1e-4) # smaller step sizes for search
               )))
}


# Zero-inflated beta-binomial model (two parameters) log-likelihood
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # param: two-element vector c(psi, lambda)
ZI.BB2.logL <- function(x, n, param){
  
  psi <- param[1]
  lambda <- param[2]
  
  require(rmutil)
  
  m <- 1 / (2 + lambda)
  s <- 2 + lambda
  
  return(sum(log((x == 0) * (1 - psi) + dbetabinom(x, n, m, s) * psi)))
}

# Fit zero-inflated beta-binomial model (two parameters) via MLE
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # psi.start: initial guess for psi
    # lambda.start: initial guess for lambda
ZI.BB2.mle <- function(x, n, psi.start, lambda.start){
  
  return(optim(par = c(psi.start, lambda.start), 
               fn = ZI.BB2.logL,
               method = "L-BFGS-B",
               x = x,
               n = n,
               lower = c(1e-6, 0),
               upper = c(1, 1e6),
               control = list(
                 fnscale = -1
               )))
}


# Zero-inflated beta-binomial model (three parameters) log-likelihood
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # param: three-element vector c(psi, alpha, beta)
ZI.BB3.logL <- function(x, n, param){
  
  psi <- param[1]
  m <- param[2] / (param[2] + param[3])
  s <- param[2] + param[3]
  
  require(rmutil)
  
  return(sum(log((x == 0) * (1 - psi) + dbetabinom(x, n, m, s) * psi)))
}

# Fit zero-inflated beta-binomial model (three parameters) via MLE
  # ARGUMENTS:
    # x: vector of number of successes (# Paranthropus specimens)
    # n: vector of number of trials (# total mammalian specimens)
    # psi.start: initial guess for psi
    # alpha.start: initial guess for alpha
    # beta.start: initial guess for beta
ZI.BB3.mle <- function(x, n, psi.start, alpha.start, beta.start){
  
  return(optim(par = c(psi.start, alpha.start, beta.start), 
               fn = ZI.BB3.logL,
               method = "L-BFGS-B",
               x = x,
               n = n,
               lower = c(1e-6, 1e-6, 1e-6),
               upper = c(1, 1e10, 1e10),
               control = list(
                 fnscale = -1
               )))
}


# function for calculating AICc (Burnham & Anderson 2002, pg 324)
  # ARGUMENTS:
    # logL: log-likelihood value
    # k: number of parameters
    # n: sample size
AICc <- function(logL, k, n){
  
  return(-2 * logL + 2 * k + 2 * k * (k + 1) / (n - k - 1))
}


# function for calculating weights from AICc scores (Burnham & Anderson 2002, pg 75)
  # ARGUMENTS:
    # AICc.res: vector of AICc calculations from different models
AICc_weights <- function(AICc.res){
  
  AICc.delta <- AICc.res - min(AICc.res)
  
  return(exp(-0.5 * AICc.delta) / sum(exp(-0.5 * AICc.delta)))
}