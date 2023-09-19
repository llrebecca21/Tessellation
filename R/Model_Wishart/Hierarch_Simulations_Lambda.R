# Multiple Single Stationary Time Series with unstructured covariance (\Lambda^{-1})
library(mvtnorm)
library(progress)
library(fda)
library(bayesplot)
library(ggplot2)
library(mcmcse)
set.seed(100)
source("R/Model_Wishart/posterior_hierarch_Lambda.R")
source("R/Model_Wishart/gradient_hierarch_Lambda.R")
source("R/Model_Wishart/he_hierarch_Lambda.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")
source("R/Model_Wishart/Sim_hierarch_Lambda.R")

# Set outer parameters for simulations
n = 1000
iter = 1000
R = 8
B = 10
J = floor((n-1) / 2)
burnin = 10

# How many simulations:
Sim = 1

#
bb_beta_big = array(data = NA, dim = c(Sim,iter,B+1,R))
Lambda_big = array(data = NA, dim = c(Sim, iter, B+1, B+1))
Theta_big = array(data = NA, dim = c(Sim, iter, B+2))
theta_true_big = array(data = NA, dim = c(Sim, R))
perio_big = array(data = NA, dim = c(Sim, J+1 , R))



pb = progress_bar$new(total = Sim)
t1 = Sys.time()
for(i in 1:Sim){
  funccall = Sim_hierarch_Lambda(n = n, iter = iter, R = R, B = B)
  # extract results
  bb_beta_big[i,,,] = funccall$bb_beta_array
  Lambda_big[i,,,] = funccall$Lambda_array
  Theta_big[i,,] = funccall$Theta
  theta_true_big[i,] = funccall$theta_vec
  perio_big[i,,] = funccall$perio
  pb$tick()
}
Sys.time() - t1

# Define omega and Psi for plotting purposes
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

theta_true_big 

# Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4))
for(r in 1:R){
  specdens = exp(Psi %*% t(bb_beta_big[1,,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-4,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  for(h in sample(ncol(specdens), 100, replace = FALSE)){
    lines(x = omega, y = log(specdens[,h]), col = rgb(0, 0, 0, 0.2))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true_big[1,r])), col = "red", lwd = 2)
  points(x = omega, y = log(perio_big[1,,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Estimate", "True"))
}


# Plot the posterior mean log spectral density against the true for a single simulation
# Grab one of the Simulation
s = sample(Sim, 1)
# get true specdens
par(mfrow = c(1,1))
true_logspec = log(arma_spec(omega = omega, theta = theta_true_big[s,r]))
plot(true_logspec, x = omega, type = "l", lwd = 2, col = "red")

# Calculate the posterior mean of the log spectral density of the estimates
specdens_mean = rowMeans(Psi %*% t(bb_beta_big[s,,,r]))
lines(specdens_mean, x = omega)

# Create for loop to plot each 
# Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4))
for(r in 1:R){
  specdens = rowMeans(Psi %*% t(bb_beta_big[s,,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-4,2), ylab = "Spectral Density", xlab = "omega",
       main = sprintf("Spectral Density\nr = %g", r))
  lines(x = omega, y = specdens)
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true_big[s,r])), col = "red", lwd = 2)
  #points(x = omega, y = log(perio_big[s,,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Post Mean", "True"), bty = "n")
}








