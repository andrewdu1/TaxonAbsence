# Function to fit zero-inflated beta-binomial model using Metropolis-Hastings algorithm

# beta distribution within beta-binomial is specified with one-parameter (lambda) and is monotonically decreasing 

# prior on lambda follows a gamma distribution since its support is all non-negative real numbers

# prior on psi is a beta distribution

# Data
  # y: vector of number of Paranthropus specimens across sites
  # n: vector of number of total mammalian specimens across sites

# Hyperparameters
  # psi.alpha: alpha parameter for beta distribution
  # psi.beta: beta parameter for beta distribution
  # lambda.shape: shape parameter for gamma distribution
  # lambda.rate: rate parameter for gamma distribution

# Initial guesses
  # psi.init: support -> real numbers from zero to one
  # lambda.init: support -> all non-negative real numbers
  # z.init: vector of latent indicator variables (0 or 1)

# Misc.
  # lambda.tune: SD tuning parameter for lambda proposal distribution
  # n.mcmc: number of MCMC iterations
  # silent: if FALSE, return progress of MCMC algorithm

ZI.BB.mcmc <- function(
    y, 
    n, 
    psi.alpha, 
    psi.beta, 
    lambda.shape, 
    lambda.rate, 
    psi.init, 
    lambda.init, 
    z.init, 
    lambda.tune, 
    n.mcmc, 
    silent = FALSE){
  
  library(VGAM) # for beta-binomial functions
  library(truncnorm) # for truncated normal functions
  
  # empty arrays to save results to
  psi.save <- lambda.save <- numeric(n.mcmc)
  z.save <- matrix(nrow = length(y), ncol = n.mcmc)
  
  # save initial values to another object to be written over
  psi <- psi.init
  lambda <- lambda.init
  z <- z.init
  
  # MCMC loop
  for(k in seq_len(n.mcmc)){
    
    # sample z (for y==0 only)
    theta <- beta(1, n[y == 0] + 1 + lambda) / beta(1, 1 + lambda)
    p.tilde <- psi * theta / (psi * theta + 1 - psi)
    z[y == 0] <- rbinom(sum(y == 0), 1, p.tilde) 
    
    # sample psi
    psi <- rbeta(1, sum(z) + psi.alpha, sum(1 - z) + psi.beta)
    
    # sample lambda (for z==1 only)
    ## sample lambda from proposal distribution
    lambda.star <- rtruncnorm(
      n = 1,
      a = 0,
      b = Inf,
      mean = lambda,
      sd = lambda.tune
    )
    ## calculate M-H ratio
    mh1 <- sum(
      dbetabinom.ab(
        x = y[z == 1], 
        size = n[z == 1],
        shape1 = 1,
        shape2 = 1 + lambda.star,
        log = TRUE)) + 
      dgamma(
        x = lambda.star,
        shape = lambda.shape,
        rate = lambda.rate,
        log = TRUE)
    mh2 <- sum(
      dbetabinom.ab(
        x = y[z == 1], 
        size = n[z == 1],
        shape1 = 1,
        shape2 = 1 + lambda,
        log = TRUE)) + 
      dgamma(
        x = lambda,
        shape = lambda.shape,
        rate = lambda.rate,
        log = TRUE)
    mh <- exp(mh1 - mh2)
    
    if(mh > runif(1)) lambda <- lambda.star
    
    # save results
    z.save[, k] <- z
    psi.save[k] <- psi
    lambda.save[k] <- lambda
    
    # keep track of progress
    if(!silent & k %% 1000 == 0) cat(k, " ") 
  }
  
  # return output
  return(list(z = z.save, psi = psi.save, lambda = lambda.save))
}