# Multiple Single Stationary Time Series from the same underlying spectra
library(mvtnorm)
library(progress)
library(fda)
library(bayesplot)
library(ggplot2)
set.seed(31)
source("R/posterior_multiple.R")
source("R/gr_multiple.R")
source("R/he_multiple.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")
source("R/mcmc_stationary.R")

# Set outer parameters for simulations
n = 1000
iter = 1000
R = 10
B = 10

# create vector of frequencies to check
test_omega = seq(0, pi, length.out = 10)
# True spectral density value at a given omega
true_spec = arma_spec(omega = test_omega, phi = 0.5, theta = 0)
# Define Psi
test_Psi = outer(X = test_omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
test_Psi[,1] = 1

# Obtaining "true" Betas
full_omega = seq(0, pi, length.out = 100000)
full_spec = arma_spec(omega = full_omega, phi = 0.5, theta = 0)
full_Psi =  outer(X = full_omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
full_Psi[,1] = 1
true_Beta = coef(lm(log(full_spec) ~ full_Psi + 0))

# How many simulations:
Sim = 500

# Create object to store coverage checks
Coverage = matrix(NA, nrow = Sim, ncol = length(test_omega))

# Create object to store simultaneous beta checks
Beta_distance = rep(NA, Sim)

pb = progress_bar$new(total = Sim)
t1 = Sys.time()
# construct for loop for simuation:
for(i in 1:Sim){
  
  pb$tick()
  
  Theta = mcmc_stationary(n = n, iter = iter, R = R, B = B)
  
  test_specdens = exp(test_Psi %*% t(Theta[ ,-(B+2)]))
  
  # Obtain the quantiles from specdens
  crit_bounds = apply(test_specdens, 1, FUN = function(x){quantile(x,c(0.025, 0.975))})
  
  # Check if true_spec is in crit_bounds
  coverage_check = ifelse(true_spec <= crit_bounds[2,] & true_spec >= crit_bounds[1,], yes = 1, no = 0)
  
  # For point-wise for spectral density
  Coverage[i,] = coverage_check
  
  # Simultaneous for Beta
  Beta = Theta[, -c(B+2)]
  # Center Beta
  Beta = Beta - tcrossprod(rep(1,nrow(Beta)), colMeans(Beta))
  # Estimate covariance matrix of Beta (sample covariance matrix)
  Beta_cov = crossprod(Beta) / nrow(Beta)
  # Calculate the weighted distance of true beta and sample beta
  x = (true_Beta - colMeans(Theta[,-(B+2)]))
  Beta_distance[i] = c(t(x) %*% solve(Beta_cov, x))
  
  
}
Sys.time() - t1


# Point-wise Coverage estimate
colMeans(Coverage)

# Simultaneous Beta Coverage
print(Beta_distance)
hist(Beta_distance)


# Simultaneous Beta distance critical value
# qchisq(p = .95, df = B + 1)

mean(Beta_distance <= qchisq(p = .95, df = B + 1))




#################
# Sim Results
#################
# Inputs: 
# Run Time : 1.573408 mins
# n = 500
# iter = 1000
# R = 1
# B = 10
# Sim = 100
# Simultaneous Coverage :
# 0.82
# Pointwise Coverage :
# 0.70 0.75 0.96 0.91 0.89 0.87 0.88 0.90 0.92 0.93


# Inputs :
# Run Time : 8.02804 mins
# n = 500
# iter = 1000
# R = 1
# B = 10
# Sim = 500
# Simulatenous Coverage :
# 0.832
# Pointwise Coverage :
# 0.722 0.774 0.912 0.914 0.878 0.908 0.918 0.918 0.896 0.914


# Inputs :
# Run Time : 9.848147 mins
# n = 500
# iter = 1000
# R = 5
# B = 10
# Sim = 500
# Simultaneous Coverage :
# 0.674
# Pointwise Coverage :
# 0.810 0.854 0.880 0.864 0.882 0.880 0.880 0.876 0.872 0.888


# Inputs :
# Run Time : 5.176307 mins
# n = 500
# iter = 1000
# R = 5
# B = 5
# Sim = 500
# Simultaneous Coverage:
# 0.61
# Pointwise Convergage :
# 0.722 0.774 0.912 0.914 0.878 0.908 0.918 0.918 0.896 0.914


# Inputs :
# Run Time : 20.54113 mins (even with computer pausing for a minute or so)
# n = 1000
# iter = 1000
# R = 10
# B = 10
# Sim = 500
# Simultaneous Coverage:
# 0.464
# Pointwise Convergage :
# 0.812 0.842 0.840 0.852 0.838 0.844 0.852 0.848 0.866 0.818

# Try different using Fisher's F-ratio
Frat = nrow(Beta) * (Beta_distance) * (nrow(Beta) - B - 1) / ((B+1) *(nrow(Beta) - 1) * (nrow(Beta) + 1))
critf = qf(0.95, B+1, nrow(Beta) - (B+1))
mean(Frat <= critf)
# 0.486






