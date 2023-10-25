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
X <- d$Paran_nisp # Paranthropus abundance
n <- X + d$NonParanMamm_nisp # Total large mammalian abundance

# Fit model using expectation-maximization algorithm
EM.res <- EM(X = X, n = n, psi.init = 0.5, lambda.init = 10)

# Get out estimated parameters
psi_hat <- EM.res$psi_hat
lambda_hat <- EM.res$lambda_hat

# Calculate posterior probabilities (prob. Paranthropus is truly absent given estimated parameters and data)
tau_hat <- 1 - tau(X, n, psi_hat, lambda_hat)

tau.df <- data.frame(site = paste(d$Site_Formation, d$Member, sep = "_"), tau_hat)

tau.df <- tau.df[tau.df$tau_hat > 0, ] # remove sites where Paranthropus is present

# save EM model fitting results
#saveRDS(EM.res, file = "Results/EM results.rds")


#############################################

# bootstrap 95% CIs for hyperparameters
set.seed(99)

n.iter <- 1000 # number of bootstrap iterations

psi_boot <- lambda_boot <- numeric(n.iter) # create empty vectors to save bootstrapped values to

for(i in seq_along(psi_boot)){
  
  index <- sample.int(length(X), size = length(X), replace = TRUE) # randomly choose index by which to resample sites
  
  X_boot <- X[index] # get resampled X
  n_boot <- n[index] # get resampled n
  
  X_boot_spec <- numeric(length(X)) # empty vector to save bootstrapped specimens to
  
  for(j in seq_along(X_boot_spec)) X_boot_spec[j] <- rbinom(1, n_boot[j], prob = X_boot[j] / n_boot[j]) # for each bootstrapped X and n, resample X (number of specimens)
  
  EM.boot <- EM(X = X_boot_spec, n = n_boot, psi.init = 0.5, lambda.init = -5) # use EM to estimate hyperparameters with bootstrapped X and n
  
  # save results to empty vectors
  psi_boot[i] <- EM.boot$psi_hat
  lambda_boot[i] <- EM.boot$lambda_hat
  
  print(i)
}

# 95% CIs
quantile(psi_boot, c(0.025, 0.975)) # 0.29-0.72
quantile(lambda_boot, c(0.025, 0.975)) # -367, -65

# save bootstrapped results
#saveRDS(list(psi_boot = psi_boot, lambda_boot = lambda_boot), file = "Results/Bootstrap CI results.rds")


#############################################

# Calculated expected frequency of Paranthropus specimens across sites using estimated model parameters
BBpdf_expect <- lapply(n, function(n1) stdReflectBBpdf(seq(0, n1), n1, lambda_hat)) # calculate probabilities for sampling number of Paranthropus specimens from zero all the way to total number of mammalian specimens at site. Do this using the beta-binomial PMF.

BBpdf_expect <- lapply(BBpdf_expect, function(x) c(x, rep(0, (max(n) + 1) - length(x)))) # add zeros to the end of each vector for NISP that exceeds the total found at each site

BBpdf_expect_sum <- Reduce("+", BBpdf_expect) # sum PMF vectors across sites

X_expect <- BBpdf_expect_sum * psi_hat # multiply summed probabilities by probability that Paranthropus is at site (psi_hat)

X_expect[1] <- X_expect[1] + (1 - psi_hat) * length(BBpdf_expect) # add (1 - psi_hat) times number of sites to the expected number of sites with zero Paranthropus. (1 - psi_hat) is the probability Paranthropus is not at the site.

# save Paranthropus expected frequency distribution results
#saveRDS(X_expect, file = "Results/Expected freq dist Paran NISP.rds")