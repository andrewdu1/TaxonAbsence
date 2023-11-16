###########################
#                         #
#     Model selection     #
#                         #
###########################

# Author: Andrew Du

# This script uses AICc to assess our model's fit (ZI.BB) against other models


# source functions
source("Code/1_R functions.R")

# Read in data
d <- read.csv("Datasets/NISP data.csv")
x <- d$Paran_nisp # Paranthropus NISP
n <- d$Paran_nisp + d$NonParanMamm_nisp # mammalian NISP

# create empty list to save model fitting results to
mle.res <- vector(
  mode = "list", 
  length = 5)

## name list elements
names(mle.res) <- c(
  "binom",
  "BB",
  "ZI.binom",
  "ZI.BB",
  "ZI.BB2"
)

# Fit models
mle.res$binom <- binom.mle(x, n, p.start = 0.5)
mle.res$BB <- BB.mle(x, n, lambda.start = 100)
mle.res$ZI.binom <- ZI.binom.mle(
  x,
  n,
  psi.start = 0.5,
  p.start = 0.001,
  psi.bounds = c(1e-6, 1),
  p.bounds = c(1e-4, 0.02)
)
mle.res$ZI.BB <- ZI.BB.mle(
  x, 
  n, 
  psi.start = 0.5, 
  lambda.start = 100
)
mle.res$ZI.BB2 <- ZI.BB2.mle(
  x,
  n,
  psi.start = 0.5,
  alpha.start = 5,
  beta.start = 5
)

# NB: ZI.BB2 gave nonsensical parameter estimates (i.e., psi_hat = 1),
# so we excluded this model from the model selection procedure, and 
# we explored its behavior via simulations
mle.res1 <- mle.res[names(mle.res) != "ZI.BB2"]

# Calculate AICc
AICc.res <- sapply(mle.res1, function(x){
  AICc(logL = x$value, k = length(x$par), n = nrow(d))
})

# Calculate AICc weights
AICc_weights(AICc.res)

# Save results
saveRDS(list(
    logL = sapply(mle.res1, function(x) x$value),
    k = sapply(mle.res1, function(x) length(x$par)),
    AICc = AICc.res, 
    AICc_w = AICc_weights(AICc.res)
    ), 
  file = "Results/Model selection results.rds")
