# Multiple Stationary Time Series simulated from perturbations with the same mean from a single time series
library(mvtnorm)
library(progress)
library(fda)
library(bayesplot)
library(ggplot2)
source("R/General_Functions/Chol_sampling.R")
source("R/General_Functions/arma_spec.R")
source("R/Data_Generation/data_generation.R")
source("R/Model_eta_br/gradient_eta_br.R")
source("R/Model_eta_br/he_eta_br.R")
source("R/Model_eta_br/posterior_eta_br.R")
# Set outer parameters for simulations
n = 1000 #time series length
iter = 1000
R = 8 #number of time series
B = 10 #number of basis coefficients
# highest little j index value for the frequencies
J = floor((n-1) / 2)
# Frequency (\omega_j): defined on [0, 2\pi)
omega = (2 * pi * (0:J)) / n

# generate data
gendata = generate_Krafty(n = n, R = R)
timeseries = gendata$matrix_timeseries
theta_true = gendata$theta_true

# Define Periodogram
# Define y_n(\omega_j) for the posterior function below
perio = (abs(mvfft(timeseries)) ^ 2 / n)
# subset perio for unique values, J = ceil((n-1) / 2) 
perio = perio[(0:J) + 1, , drop=FALSE]


##########################
# Set Hyper-Parameters
##########################
# Values that are fixed and set by user
# Create matrix of the basis functions
# fix fourier frequencies
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1
# Define degrees of freedom for half-t prior: \pi(\tau^2 | \lambda)
# Cauchy nu = 1
nu = 3
# Define eta as the other scale parameter for half-t prior: \pi(\lambda)
# etasq = 1 gives standard Cauchy; higher eta gives wider Cauchy
etasq = 1
# Define D's main diagonal : 
# D is a measure of prior variance for \beta_1 through \beta_K
# Rebecca's D
D = 1 / (4 * pi * (1:B)^2)

# Calculate ybar(omega_j)
y_bar = rowMeans(perio)

########################
# Initialize Parameters
########################
# set tau^2 value
tausquared = 50
# Define lambda for the half-t prior: \pi(\tau^2 | \lambda) and \pi(\lambda)
lambda = 1
# The new D matrix that houses the prior variance of \beta^* 
Sigma = c(1, D * tausquared)
# Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
# ncol: number of parameters (beta^*, tau^2)
# dim : (iter) x (B + 2)
Theta = matrix(NA, nrow = iter, ncol = B + 2)
# Using J amount of data for periodogram, can initialize beta this way:
betavalues = solve(crossprod(Psi), crossprod(Psi, log(y_bar)))
# Initialize bb_beta at the mean for the prior of bb_beta
bb_beta = tcrossprod(betavalues, rep(1,R))
# Specify Sum of X for the posterior function later
# Specify Sum of X for the posterior function later
# 1^T_n X part in the paper: identical to colSums but is a faster calculation
sumPsi = crossprod(Psi, rep(1,J+1)) 
# Initialize first row of Theta
Theta[1,] = c(betavalues, tausquared)
# Create Array to store values of bb_beta
bb_beta_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
# initialize first array with bb_beta value
bb_beta_array[1,,] = bb_beta
# initialize starting value of eta_br
eta_br = matrix(1, nrow = B + 1 , ncol = R)
# Create Array to store eta_br values
eta_br_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
# Initialize eta_br values
eta_br_array[1,,] = eta_br

#####################
# MCMC Algorithm
#####################
# need to update blackboard beta, eta_b^r, tau^2, lambda 

pb = progress_bar$new(total = iter - 1)
t1 = Sys.time()
for (g in 2:iter) {
  pb$tick()
  #g = 2
  ###############################################
  # Sample \tau^2 and \lambda jointly with Gibbs
  ###############################################
  # Update \lambda and \tau with Gibbs Sampler
  lambda = 1/rgamma(1, (nu+1)/2, nu/tausquared + etasq)
  tausquared = 1/rgamma(1, (nu + B * R)/2, sum(rowSums(bb_beta^2 * eta_br)[-1] / D)/2 + nu/lambda)
  # Update Theta matrix with new tau squared value
  Theta[g,B+2] = tausquared
  # Update Sigma with new tau^2 value
  Sigma = c(1, D * tausquared)
  #######################################
  # Sample \eta_b^r using slicing method
  #######################################
  rate = c(bb_beta * bb_beta / (2 * Sigma))
  p = R * (B + 1)
  u = runif(p)/(1 + c(eta_br))
  q = runif(p)*(pgamma(rate/u^2, shape = 1))
  eta_br = qgamma(q, shape = 1)/rate
  # reconstruct shape
  eta_br = matrix(eta_br, nrow = B + 1 , ncol = R)
  # put eta_br into array for storage
  eta_br_array[g,,] = eta_br

  ######################
  # bb_beta update :MH
  ######################
  # Update \mathbb{B} with Metropolis Hastings Sampler
  for(r in 1:R){
    #r = 1
    br = bb_beta[,r]
    eta_r = eta_br[,r]
    # Maximum A Posteriori (MAP) estimate : finds the \beta that gives us the mode of the conditional posterior of \beta conditioned on y
    map = optim(par = br, fn = posterior_eta_br, gr = gradient_eta_br, method ="BFGS", control = list(fnscale = -1),
                 Psi = Psi, sumPsi = sumPsi, y = perio[,r], Sigma = Sigma, eta_r = eta_r)$par
    # Call the hessian function
    norm_precision = he_eta_br(br = map, Psi = Psi, y = perio[,r], Sigma = Sigma, eta_r = eta_r) * -1
    # Calculate the \beta^* proposal, using Cholesky Sampling
    betaprop = Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
    # Calculate acceptance ratio
    prop_ratio = min(1, exp(posterior_eta_br(br = betaprop,  Psi = Psi, sumPsi = sumPsi, y = perio[,r], eta_r = eta_r, Sigma = Sigma) -
                               posterior_eta_br(br = br, Psi = Psi, sumPsi = sumPsi, y = perio[,r], eta_r = eta_r, Sigma = Sigma)))
    # Create acceptance decision
    accept <- runif(1)
    if(accept < prop_ratio){
      # Accept betaprop as new beta^(r)
      bb_beta_array[g, ,r] <- betaprop
    }else{
      # Reject betaprop as new beta^(r)
      bb_beta_array[g, ,r] <- br
    }
  }
  
}
Sys.time() - t1

#######################
# Plots and Diagnostics
#######################
# Need to update this code for the hierarchical set-up
# Remove burn-in
burnin = 10
Theta = Theta[-(1:burnin),]
bb_beta_array = bb_beta_array[-(1:burnin),,]
eta_br_array = eta_br_array[-(1:burnin),,]
r = 8


# Create spectral densities for the first series
specdens = exp(Psi %*% t(bb_beta_array[,,r]))
par(mfrow = c(1, 1))
plot(x =c(), y=c(), xlim = c(0,3), ylim = range(log(specdens)), ylab = "Spectral Density", xlab = "omega",
     main = "Spectral Density Estimates \nwith True Spectral Density")
for(h in sample(ncol(specdens), 100, replace = FALSE)){
  lines(x = omega, y = log(specdens[,h]), col = rgb(0, 0, 0, 0.2))
}
lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "red", lwd = 2)
#abline(v = pi/4)
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Estimate", "True"))


theta_true

Lambda_array[900,,]



# Plot the Spectral Density Estimates
#pdf(file = "Spectral_Density_Estimates.pdf",
#    width = 10,
#    height = 5,)
specdens = exp(Psi %*% t(Theta[ ,-(B+2)]))
par(mfrow = c(1, 1))
plot(x =c(), y=c(), xlim = c(0,3), ylim = range(log(specdens)), ylab = "Spectral Density", xlab = "omega",
     main = "Spectral Density Estimates \nwith True Spectral Density")
for(h in sample(ncol(specdens), 1000, replace = FALSE)){
  lines(x = omega, y = log(specdens[,h]), col = rgb(0, 0, 0, 0.2))
}
lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "red", lwd = 2)
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Estimate", "True"))
#dev.off()

dim(specdens)
# n X (numiter - burn)

# Create Data frame to store the lower bound and upper bound and mean
summary_stats <- data.frame("lower" = apply(specdens, 1, FUN = function(x){quantile(x, .025)}), "mean" = rowMeans(specdens),
                            "upper" = apply(specdens, 1, FUN = function(x){quantile(x, 0.975)}))



















