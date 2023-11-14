# Read in simulation results
BB.sim.res <- readRDS("Results/BB sim results.rds")

# get out simulation parameters
psi.sim <- as.numeric(dimnames(BB.sim.res)[[1]])
lambda.sim <- as.numeric(dimnames(BB.sim.res)[[2]])
n.iter <- dim(BB.sim.res)[3]
param <- dimnames(BB.sim.res)[[4]]

# number of different combinations of parameter estimates
BB.nrow <- prod(dim(BB.sim.res))

# create empty data frame
BB.sim <- data.frame(
  estimate = rep(NA, BB.nrow),
  parameter = rep(NA, BB.nrow),
  psi = rep(NA, BB.nrow),
  lambda = rep(NA, BB.nrow)
)

# fill out parameter column
BB.sim$parameter <- rep(param, each = n.iter)

# get all combinations of psi & lambda
param.comb <- expand.grid(psi.sim, lambda.sim)

# fill out psi & lambda column
BB.sim[, c("psi", "lambda")] <- 
  param.comb[
    rep(
      seq_len(nrow(param.comb)), 
      times = n.iter
      ), 
    ]

# fill out estimates column
for(psi_i in seq_along(psi.sim)){
  for(lambda_i in seq_along(lambda.sim)){
    for(j in seq_along(param)){
      
      BB.sim$estimate[
        BB.sim$parameter == param[j] &
        BB.sim$psi == psi.sim[psi_i] & 
        BB.sim$lambda == lambda.sim[lambda_i]
      ] <- BB.sim.res[psi_i, lambda_i, , j]
    }
  }
}

# load ggplot
library(ggplot2)

#pdf("Figures/ZI-BB simulation results.pdf", height = 9)

# psi_hat plot
BB.sim.psi <- BB.sim[BB.sim$parameter == "psi_hat", ]

## change facet labels
psi_labs <- paste("True psi =", psi.sim)
names(psi_labs) <- psi.sim

## create plot
ggplot(BB.sim.psi, aes(x = factor(lambda), y = estimate)) + 
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
    data = BB.sim.psi, 
    aes(yintercept = psi), 
    linetype = "dashed") + 
  theme_bw() +
  labs(
    x = expression("True lambda (" * lambda * ")"),
    y = expression("Estimated psi (" * hat(psi) * ")")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))
  

# lambda_hat plot
BB.sim.lambda <- BB.sim[BB.sim$parameter == "lambda_hat", ]

## change facet labels
lambda_labs <- paste("True lambda =", lambda.sim)
names(lambda_labs) <- lambda.sim

## create plot
ggplot(BB.sim.lambda, aes(x = factor(psi), y = estimate)) + 
  geom_violin(trim = FALSE, fill = "lightskyblue1") +
  stat_summary(
    fun.y = mean, 
    geom = "point", 
    shape = 23, 
    size = 2) +
  facet_wrap(
    ~lambda, 
    ncol = 2, 
    labeller = labeller(lambda = lambda_labs),
    scales = "free_y") +
  geom_hline(
    data = BB.sim.lambda, 
    aes(yintercept = lambda), 
    linetype = "dashed") + 
  theme_bw() +
  labs(
    x = expression("True psi (" * psi * ")"),
    y = expression("Estimated lambda (" * hat(lambda) * ")")) +
  scale_y_log10() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12))

#dev.off()