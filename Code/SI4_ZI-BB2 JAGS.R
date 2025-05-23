# Function to fit zero-inflated beta-binomial model using JAGS

# beta distribution within beta-binomial is specified with one parameter (lambda) and is monotonically decreasing 

# hyperprior on lambda follows a gamma distribution since its support is all non-negative real numbers

# hyperprior on psi is a beta distribution

# psi: marginal probability of presence
# lambda: controls shape of beta distribution in beta-binomial
# p: sampling probability
# y: number of Paranthropus specimens
# n: number of total mammalian specimens

model{
  # hyperpriors
  psi ~ dbeta(1, 1)
  lambda ~ dgamma(0.001, 0.001)
  
  for(i in 1:length(y)){
    # prior
    p[i] ~ dbeta(1, 1 + lambda)
    # likelihood
    z[i] ~ dbern(psi)
    y[i] ~ dbin(z[i] * p[i], n[i])
  }
}