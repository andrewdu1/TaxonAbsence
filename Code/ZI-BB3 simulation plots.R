# Read in simulation results
BB3.sim.res <- readRDS("Results/BB3 sim results.rds")

# get out simulation parameters
psi.sim <- as.numeric(dimnames(BB3.sim.res)[[1]])
alpha.sim <- as.numeric(dimnames(BB3.sim.res)[[2]])
beta.sim <- as.numeric(dimnames(BB3.sim.res)[[3]])
n.iter <- dim(BB3.sim.res)[4]
param <- dimnames(BB3.sim.res)[[5]]

# number of different combinations of parameter estimates
BB3.nrow <- prod(dim(BB3.sim.res))

# create empty data frame
BB3.sim <- data.frame(
  estimate = rep(NA, BB3.nrow),
  parameter = rep(NA, BB3.nrow),
  psi = rep(NA, BB3.nrow),
  A = rep(NA, BB3.nrow),
  B = rep(NA, BB3.nrow)
)

# fill out parameter column
BB3.sim$parameter <- rep(param, each = BB3.nrow / length(param))

# get all combinations of psi, alpha, & beta
param.comb <- expand.grid(
  psi.sim, 
  alpha.sim, 
  beta.sim)

# fill out psi, alpha, & beta column
BB3.sim[, c("psi", "A", "B")] <- 
  param.comb[
    rep(
      seq_len(nrow(param.comb)), 
      times = n.iter
    ), 
  ]

# fill out estimates column
for(psi_i in seq_along(psi.sim)){
  for(alpha_i in seq_along(alpha.sim)){
    for(beta_i in seq_along(beta.sim)){
      for(j in seq_along(param)){
        
        BB3.sim$estimate[
          BB3.sim$parameter == param[j] &
          BB3.sim$psi == psi.sim[psi_i] & 
          BB3.sim$A == alpha.sim[alpha_i] &
          BB3.sim$B == beta.sim[beta_i]
        ] <- BB3.sim.res[psi_i, alpha_i, beta_i, , j]
      }
    }
  }
}

# load ggplot
library(ggplot2)

#pdf("Figures/ZI-BB3 simulation results.pdf", height = 9)

# psi_hat plot
BB3.sim.psi <- BB3.sim[BB3.sim$parameter == "psi_hat", ]

## change facet labels
psi_labs <- paste("True psi =", psi.sim)
names(psi_labs) <- psi.sim

## create plot
### psi ~ alpha
ggplot(BB3.sim.psi, aes(x = factor(A), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "lightpink1") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~psi, 
    ncol = 2, 
    labeller = labeller(psi = psi_labs)) +
  geom_hline(
    data = BB3.sim.psi, 
    aes(yintercept = psi), 
    linetype = "dashed") + 
  theme_bw() +
  labs(
    x = expression("True alpha (" * alpha * ")"),
    y = expression("Estimated psi (" * hat(psi) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))

### psi ~ beta
ggplot(BB3.sim.psi, aes(x = factor(B), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "lightpink3") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~psi, 
    ncol = 2, 
    labeller = labeller(psi = psi_labs)) +
  geom_hline(
    data = BB3.sim.psi, 
    aes(yintercept = psi), 
    linetype = "dashed") + 
  theme_bw() +
  labs(
    x = expression("True beta (" * beta * ")"),
    y = expression("Estimated psi (" * hat(psi) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))


# alpha_hat plot
BB3.sim.alpha <- BB3.sim[BB3.sim$parameter == "alpha_hat", ]

## change facet labels
alpha_labs <- paste("True alpha =", alpha.sim)
names(alpha_labs) <- alpha.sim

## create plot
### alpha ~ psi
ggplot(BB3.sim.alpha, aes(x = factor(psi), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "seagreen1") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~A, 
    ncol = 2, 
    labeller = labeller(A = alpha_labs)) +
  geom_hline(
    data = BB3.sim.alpha, 
    aes(yintercept = A), 
    linetype = "dashed") + 
  scale_y_log10() + 
  theme_bw() +
  labs(
    x = expression("True psi (" * psi * ")"),
    y = expression("Estimated alpha (" * hat(alpha) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))

### alpha ~ beta
ggplot(BB3.sim.alpha, aes(x = factor(B), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "seagreen3") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~A, 
    ncol = 2, 
    labeller = labeller(A = alpha_labs)) +
  geom_hline(
    data = BB3.sim.alpha, 
    aes(yintercept = A), 
    linetype = "dashed") + 
  scale_y_log10() +
  theme_bw() +
  labs(
    x = expression("True psi (" * psi * ")"),
    y = expression("Estimated beta (" * hat(beta) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))


# beta_hat plot
BB3.sim.beta <- BB3.sim[BB3.sim$parameter == "beta_hat", ]

## change facet labels
beta_labs <- paste("True beta =", beta.sim)
names(beta_labs) <- beta.sim

## create plot
### beta ~ psi
ggplot(BB3.sim.beta, aes(x = factor(psi), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "mediumpurple1") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~B, 
    ncol = 2, 
    labeller = labeller(B = beta_labs),
    scales = "free_y") +
  geom_hline(
    data = BB3.sim.beta, 
    aes(yintercept = B), 
    linetype = "dashed") + 
  scale_y_log10() +
  theme_bw() +
  labs(
    x = expression("True psi (" * psi * ")"),
    y = expression("Estimated beta (" * hat(beta) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))

### beta ~ alpha
ggplot(BB3.sim.beta, aes(x = factor(A), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "mediumpurple3") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~B, 
    ncol = 2, 
    labeller = labeller(B = beta_labs),
    scales = "free_y") +
  geom_hline(
    data = BB3.sim.beta, 
    aes(yintercept = B), 
    linetype = "dashed") + 
  scale_y_log10() +
  theme_bw() +
  labs(
    x = expression("True alpha (" * alpha * ")"),
    y = expression("Estimated beta (" * hat(beta) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))

#dev.off()