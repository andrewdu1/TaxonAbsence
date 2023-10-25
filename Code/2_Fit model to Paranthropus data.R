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
mle.res <- ZI.BB.mle(
  psi.start = 0.5,
  lambda.start = 100,
  x = x,
  n = n,
  hessian = TRUE)

# Get out estimated parameters
psi_hat <- mle.res$param["psi"]
lambda_hat <- mle.res$param["lambda"]

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

post_prob.res <- post_prob(d1$n, psi_hat, lambda_hat)

post_prob.res1 <- data.frame(
  site = d1$site, 
  post_prob = round(post_prob.res, 3))


#############################################

# bootstrap 95% CIs for hyperparameters
set.seed(99)

n.iter <- 1000 # number of bootstrap iterations

psi_boot <- lambda_boot <- numeric(n.iter) # create empty vectors to save bootstrapped values to

for(i in seq_along(psi_boot)){
  
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
  
  mle.boot <- ZI.BB.mle(
    psi.start = 0.5,
    lambda.start = 100,
    x = x_boot_spec,
    n = n_boot,
    hessian = FALSE
  ) # use ML to estimate hyperparameters with bootstrapped x and n
  
  # save results to empty vectors
  psi_boot[i] <- mle.boot["psi"]
  lambda_boot[i] <- mle.boot["lambda"]
  
  print(i)
}

# 95% CIs
quantile(psi_boot, c(0.025, 0.975)) # 0.28-0.67
quantile(lambda_boot, c(0.025, 0.975)) # 64, 361


#############################################

# Calculated expected frequency of Paranthropus specimens across sites using estimated model parameters
BBpdf_expect <- lapply(n, function(n1) stdReflectBBpdf(seq(0, n1), n1, lambda_hat)) # calculate probabilities for sampling number of Paranthropus specimens from zero all the way to total number of mammalian specimens at site. Do this using the beta-binomial PMF.

BBpdf_expect <- lapply(BBpdf_expect, function(x) c(x, rep(0, (max(n) + 1) - length(x)))) # add zeros to the end of each vector for NISP that exceeds the total found at each site

BBpdf_expect_sum <- Reduce("+", BBpdf_expect) # sum PMF vectors across sites

x_expect <- BBpdf_expect_sum * psi_hat # multiply summed probabilities by probability that Paranthropus is at site (psi_hat)

x_expect[1] <- x_expect[1] + (1 - psi_hat) * length(BBpdf_expect) # add (1 - psi_hat) times number of sites to the expected number of sites with zero Paranthropus. (1 - psi_hat) is the probability Paranthropus is not at the site.

# save Paranthropus expected frequency distribution results
#saveRDS(x_expect, file = "Results/Expected freq dist Paran NISP.rds")