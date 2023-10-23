###########################
#                         #
#     Model selection     #
#                         #
###########################

# Author: Andrew Du & Eric Friedlander

# This script uses AICc to assess our model's
# fit against simpler models


# Read in data
d <- read.csv("Datasets/NISP data.csv")
x <- d$Paran_nisp # Paranthropus NISP
n <- d$Paran_nisp + d$NonParanMamm_nisp # mammalian NISP


# create empty list to save model results to
mle.res <- vector(
  mode = "list", 
  length = 4)

## name list elements
names(mle.res) <- c(
  "binom",
  "BB",
  "ZI.binom",
  "ZI.BB"
)


# Create & fit models

## Binomial model
binom.logL <- function(x, n, p){
  
  return(sum(dbinom(x, n, p, log = TRUE)))
}

### use optim() to numerically estimate parameter
binom.fit <- optim(par = 0.5, 
             fn = binom.logL,
             method = "L-BFGS-B",
             x = x,
             n = n,
             lower = 1e-6,
             upper = 1 - 1e-6,
             control = list(
               fnscale = -1
             ))

### save logL & k (# parameters)
mle.res$binom <- c(
  binom.fit$value,
  length(binom.fit$par)
  )


## Beta-binomial model
BB.logL <- function(x, n, lambda){
  
  require(rmutil)
  
  m <- 1 / (2 + lambda)
  s <- 2 + lambda

  return(sum(dbetabinom(x, n, m, s, log = TRUE)))
}

### use optim() to numerically estimate parameter
BB.fit <- optim(par = 10, 
                fn = BB.logL,
                method = "L-BFGS-B",
                x = x,
                n = n,
                lower = 0,
                upper = Inf,
                control = list(
                fnscale = -1
                ))

### save logL & k (# parameters)
mle.res$BB <- c(
  BB.fit$value,
  length(BB.fit$par)
)


## Zero-inflated binomial model
  # Argument param is a vector with two elements
    # First element is psi
    # Second element is p
ZI.binom.logL <- function(x, n, param){
  
  psi <- param[1]
  p <- param[2]
  
  return(sum(log((x == 0) * (1 - psi) + dbinom(x, n, p) * psi)))
}

### use optim() to numerically estimate parameters
ZI.binom.fit <- optim(par = c(0.5, 0.001), 
                      fn = ZI.binom.logL,
                      method = "L-BFGS-B",
                      x = x,
                      n = n,
                      lower = c(1e-6, 1e-4),
                      upper = c(1, 0.02),
                      control = list(
                        fnscale = -1
                      ))

### save logL & k (# parameters)
mle.res$ZI.binom <- c(
  ZI.binom.fit$value,
  length(ZI.binom.fit$par)
)


## Zero-inflated beta-binomial model
  # Argument param is a vector with two elements
    # First element is psi
    # Second element is lambda
ZI.BB.logL <- function(x, n, param){
  
  psi <- param[1]
  lambda <- param[2]
  
  require(rmutil)
  
  m <- 1 / (2 + lambda)
  s <- 2 + lambda
  
  return(sum(log((x == 0) * (1 - psi) + dbetabinom(x, n, m, s) * psi)))
}

### use optim() to numerically estimate parameters
ZI.BB.fit <- optim(par = c(0.5, 100), 
                   fn = ZI.BB.logL,
                   method = "L-BFGS-B",
                   x = x,
                   n = n,
                   lower = c(1e-6, 0),
                   upper = c(1, Inf),
                   control = list(
                     fnscale = -1
                   ))

### save logL & k (# parameters)
mle.res$ZI.BB <- c(
  ZI.BB.fit$value,
  length(ZI.BB.fit$par)
)


# Write functions for AICc & AICc weights
## AICc (Burnham & Anderson 2002, pg 324)
AICc <- function(logL, k, n){
  
  return(-2 * logL + 2 * k + 2 * k * (k + 1) / (n - k - 1))
}

## AICc weights
## w_i (Burnham & Anderson 2002, pg 75)
AICc_weights <- function(AICc.res){
  
  AICc.delta <- AICc.res - min(AICc.res)
  
  return(exp(-0.5 * AICc.delta) / sum(exp(-0.5 * AICc.delta)))
}


# run functions on model results
AICc.res <- sapply(mle.res, function(x){
  
  AICc(
    logL = x[1], 
    k = x[2], 
    n = nrow(d))
})

AICc_weights(AICc.res)