##########################################
#                                        #
#     Fit model to Paranthropus data     #
#                                        #
##########################################

# Author: Andrew Du & Eric Friedlander


# Source R functions
source("Code/1_R functions.R")

# Read in data
d <- read.csv("Datasets/NISP data.csv", header = TRUE)

# Define objects
x <- d$Paran_nisp # Paranthropus abundance
n <- x + d$NonParanMamm_nisp # Total large mammalian abundance

# Fit model using MLE
ZI.BB.fit <- ZI.BB.mle(
  psi.start = 0.5,
  lambda.start = 100,
  x = x,
  n = n
  )

# Get out estimated parameters
psi_hat <- ZI.BB.fit$par[1]
lambda_hat <- ZI.BB.fit$par[2]

# Calculate posterior probabilities 
# (prob. Paranthropus is truly absent given estimated parameters and data)
d1 <- data.frame(
  site = paste(
    d$Site_Formation, 
    d$Member, 
    sep = "_"), 
  x, 
  n)

d1 <- d1[d1$x == 0, ] # remove sites where Paranthropus is present

post_prob.res <- post_prob(
  n = d1$n, 
  psi = psi_hat, 
  lambda = lambda_hat)

post_prob.res1 <- data.frame(
  site = d1$site, 
  post_prob = round(post_prob.res, 3))


#############################################

# bootstrap 95% CIs for hyperparameters & posterior probabilities
set.seed(99) # for repeatability purposes

# bootstrap the data
boot.res <- boot_data( 
  x = x,
  n = n,
  n.iter = 1000,
  silent = FALSE
)

# estimate hyperparameters for each bootstrapped iteration
mle.boot.res <- sapply(boot.res, function(i){ 
  
  mle.res <- ZI.BB.mle(
    x = i$x,
    n = i$n,
    psi.start = 0.5,
    lambda.start = 100
  )
  
  param.hat <- mle.res$par
  names(param.hat) <- c("psi", "lambda")
  return(param.hat)
})

# calculate posterior probability for each bootstrapped iteration
post_prob.boot.res <- apply(mle.boot.res, 2, function(i){ 
  
  post_prob(d1$n, i["psi"], i["lambda"])
})

# 95% CIs
quantile(mle.boot.res["psi", ], c(0.025, 0.975)) # psi: 0.28-0.67
quantile(mle.boot.res["lambda", ], c(0.025, 0.975)) # lambda: 64-361
apply(post_prob.boot.res, 1, function(x) quantile(x, c(0.025, 0.975)))


#############################################

# Calculate expected frequency distribution of Paranthropus specimens across sites using estimated model parameters
BBpmf_expect <- array( # empty array to save results to
  data = 0,
  dim = c(nrow(d), max(n) + 1) # + 1 includes zero NISP
)

for(i in seq_len(nrow(BBpmf_expect))){
  
  BBpmf_expect[i, seq(1, n[i] + 1)] <- BBpmf(
    x = seq(0, n[i]), 
    n = n[i], 
    lambda = lambda_hat)
} # calculate probabilities for sampling number of Paranthropus specimens from zero all the way to total number of mammalian specimens at each site. 

BBpmf_expect_sum <- colSums(BBpmf_expect) # sum PMF vectors across sites

x_expect <- BBpmf_expect_sum * psi_hat # multiply summed probabilities by probability that Paranthropus is at site (psi_hat)

x_expect[1] <- x_expect[1] + (1 - psi_hat) * nrow(d) # add (1 - psi_hat) times number of sites to the expected number of sites with zero Paranthropus.


#############################################

# Save results
res <- list(
  estim.param = c(psi_hat, lambda_hat),
  post_prob.res = post_prob.res1,
  mle.boot = mle.boot.res,
  post_prob.boot = post_prob.boot.res,
  x_expect = x_expect
)

#saveRDS(res, "Results/Main model results.rds")
