##########################
#                        #
#     Create figures     #
#                        #
##########################

# Author: Andrew Du
# Date: 1-12-22


# Read in data & results
d <- read.csv("Datasets/NISP data.csv", header = TRUE) # raw data

EM.res <- readRDS("Results/EM results.rds") # EM results
X_expect <- readRDS("Results/Expected freq dist Paran NISP.rds") # expected Paranthropus NISP frequency distribution

# Source functions
source("Code/1_R functions.R")

# Define objects
X <- d$Paran.nisp # Paranthropus abundance
n <- X + d$Non.paran.mamm.nisp # Total large mammalian abundance

psi_hat <- EM.res$psi_hat # estimated psi parameter
lambda_hat <- EM.res$lambda_hat # estimated lambda parameter


### Fig. 2: Histogram of number of Paranthropus and large mammal specimens

#pdf("Figures/Fig 2.pdf", height = 10)

par(mfrow = c(3, 1), mar = c(4.1, 4.5, 3, 2) + 0.1)

hist(X, col = "gray", xlab = expression(italic(Paranthropus) ~ "NISP"), ylab = "Number of sites", main = "", cex.axis = 1.75, cex.lab = 2, breaks = 10)

mtext("A", at = 0, cex = 2.25)

hist(n, col = "gray", xlab = "Total large mammalian NISP", ylab = "Number of sites", main = "", cex.axis = 1.75, cex.lab = 2, breaks = 10)

mtext("B", at = 1, cex = 2.25)

hist(X / n, col = "gray", xlab = expression(italic(Paranthropus) ~ "NISP / mammalian NISP"), ylab = "Number of sites", main = "", cex.axis = 1.75, cex.lab = 2, breaks = 10)

mtext("C", at = 0, cex = 2.25)

#dev.off()


### Fig. 3: Probability density function of Paranthropus sampling probability (using lambda parameter)
x <- seq(0, 1, length.out = 100000)

#pdf("Figures/Fig 3.pdf", height = 10)

par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 2) + 0.1)

plot(x, dbeta(x, 1, 1 - lambda_hat), type = "l", lwd = 3, xlab = expression(italic(Paranthropus) * " sampling probability (" * italic(p[i]) * ")"), ylab = "Density", cex.axis = 1.5, cex.lab = 1.5, log = "x", xaxt = "n")

axis(1, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0), labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1"), cex.axis = 1.5)

text(0.04, 130, bquote(hat(lambda) == .(round(lambda_hat))), cex = 2, pos = 4)

mtext("A", at = 0.00001, cex = 2.5)

## Expected vs. observed number of Paranthropus specimens across sites
plot(table(X), xlab = expression(italic(Paranthropus) ~ "NISP"), ylab = "Number of sites", cex.lab = 1.5, cex.axis = 1.5, lwd = 3, xaxt = "n", xlim = c(0, 34), type = "n")

points(seq(0, 34), X_expect[seq(1, 35)], pch = 21, bg = "gray85", cex = 1.5)
points(table(X), lwd = 3)

axis(1, at = seq(0, 34, 2), labels = seq(0, 34, 2), cex.axis = 1.5)

legend("topright", legend = c("Observed", "Expected"), lty = c(1, NA), pch = c(NA, 21), pt.bg = c(NA, "gray85"), pt.cex = c(NA, 2), pt.lwd = c(NA, 1), bty = "n", cex = 1.5, lwd = 3)

mtext("B", at = 0, cex = 2.5)

#dev.off()


### Fig. 4: Probability Paranthropus absence curve
n1 <- seq(0, 20000)

tau0 <- tau(X = rep(0, length(n1)), n = n1, psi = psi_hat, lambda = lambda_hat) # get taus for sites where there's no Paranthropus and with increasing number of mammalian specimens

# get out NISP corresponding to probabilities of 0.5, 0.75, 0.9, 0.95, and 0.99
nisp_0.5 <- which(abs((1 - tau0) - 0.5) == min(abs((1 - tau0) - 0.5)))

nisp_0.75 <- which(abs((1 - tau0) - 0.75) == min(abs((1 - tau0) - 0.75)))

nisp_0.9 <- which(abs((1 - tau0) - 0.9) == min(abs((1 - tau0) - 0.9)))

nisp_0.95 <- which(abs((1 - tau0) - 0.95) == min(abs((1 - tau0) - 0.95)))

nisp_0.99 <- which(abs((1 - tau0) - 0.99) == min(abs((1 - tau0) - 0.99)))


#pdf("Figures/Fig 4.pdf", height = 5)

par(mar = c(5, 6.5, 4, 2) + 0.1, fig = c(0, 1, 0, 1))

plot(n1, 1 - tau0, xlab = "Number of mammalian specimens (NISP)", ylab = expression(atop("Posterior probability " * italic(Paranthropus), " was truly absent from site")), type = "l", lwd = 3, log = "x", cex.axis = 1.5, cex.lab = 1.5)

abline(v = n1[nisp_0.5], lty = 2, lwd = 1.5)
text(n1[nisp_0.5] - 6, 0.62, paste("0.5 =", n1[nisp_0.5], "NISP"), srt = 90, cex = 1)

abline(v = n1[nisp_0.75], lty = 2, lwd = 1.5)
text(n1[nisp_0.75] - 30, 0.55, paste("0.75 =", n1[nisp_0.75], "NISP"), srt = 90, cex = 1, pos = 3)

abline(v = n1[nisp_0.9], lty = 2, lwd = 1.5)
text(n1[nisp_0.9] - 100, 0.55, paste("0.9 =", n1[nisp_0.9], "NISP"), srt = 90, cex = 1, pos = 3)

abline(v = n1[nisp_0.95], lty = 2, lwd = 1.5)
text(n1[nisp_0.95] - 200, 0.55, paste("0.95 =", n1[nisp_0.95], "NISP"), srt = 90, cex = 1, pos = 3)

abline(v = n1[nisp_0.99], lty = 2, lwd = 1.5)
text(n1[nisp_0.99] - 700, 0.55, paste("0.99 =", n1[nisp_0.99], "NISP"), srt = 90, cex = 1, pos = 3)

# plot the inset (curve on arithmetic axes)
par(fig = c(0.06, 0.5, 0.4, 0.975), new = TRUE) 

plot(n1, 1 - tau0, xlab = "", ylab = "", type = "l", lwd = 2, cex.axis = 0.75)

#dev.off()