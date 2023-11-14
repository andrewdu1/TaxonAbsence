# Read in simulation results
BB2.sim.res <- readRDS("Results/BB2 sim results.rds")

# get out simulation parameters
psi.sim <- as.numeric(dimnames(BB2.sim.res)[[1]])
alpha.sim <- as.numeric(dimnames(BB2.sim.res)[[2]])
beta.sim <- as.numeric(dimnames(BB2.sim.res)[[3]])
n.iter <- dim(BB2.sim.res)[4]
param <- dimnames(BB2.sim.res)[[5]]

# number of different combinations of parameter estimates
BB2.nrow <- prod(dim(BB2.sim.res))

# create empty data frame
BB2.sim <- data.frame(
  estimate = rep(NA, BB2.nrow),
  parameter = rep(NA, BB2.nrow),
  psi = rep(NA, BB2.nrow),
  A = rep(NA, BB2.nrow),
  B = rep(NA, BB2.nrow)
)

# fill out parameter column
BB2.sim$parameter <- rep(param, each = BB2.nrow / 3)

# get all combinations of psi, alpha, & beta
param.comb <- expand.grid(
  psi.sim, 
  alpha.sim, 
  beta.sim)

# fill out psi, alpha, & beta column
BB2.sim[, c("psi", "A", "B")] <- 
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
        
        BB2.sim$estimate[
          BB2.sim$parameter == param[j] &
          BB2.sim$psi == psi.sim[psi_i] & 
          BB2.sim$A == alpha.sim[alpha_i] &
          BB2.sim$B == beta.sim[beta_i]
        ] <- BB2.sim.res[psi_i, alpha_i, beta_i, , j]
      }
    }
  }
}

# load ggplot
library(ggplot2)

#pdf("Figures/ZI-BB2 simulation results.pdf", height = 9)

# psi_hat plot
BB2.sim.psi <- BB2.sim[BB2.sim$parameter == "psi_hat", ]

## change facet labels
psi_labs <- paste("True psi =", psi.sim)
names(psi_labs) <- psi.sim

## create plot
### psi ~ alpha
ggplot(BB2.sim.psi, aes(x = factor(A), y = estimate)) + 
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
    data = BB2.sim.psi, 
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
ggplot(BB2.sim.psi, aes(x = factor(B), y = estimate)) + 
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
    data = BB2.sim.psi, 
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
BB2.sim.alpha <- BB2.sim[BB2.sim$parameter == "alpha_hat", ]

## change facet labels
alpha_labs <- paste("True alpha =", alpha.sim)
names(alpha_labs) <- alpha.sim

## create plot
### alpha ~ psi
ggplot(BB2.sim.alpha, aes(x = factor(psi), y = estimate)) + 
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
    data = BB2.sim.alpha, 
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
ggplot(BB2.sim.alpha, aes(x = factor(B), y = estimate)) + 
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
    data = BB2.sim.alpha, 
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
BB2.sim.beta <- BB2.sim[BB2.sim$parameter == "beta_hat", ]

## change facet labels
beta_labs <- paste("True beta =", beta.sim)
names(beta_labs) <- beta.sim

## create plot
### beta ~ psi
ggplot(BB2.sim.beta, aes(x = factor(psi), y = estimate)) + 
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
    data = BB2.sim.beta, 
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
ggplot(BB2.sim.beta, aes(x = factor(A), y = estimate)) + 
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
    data = BB2.sim.beta, 
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