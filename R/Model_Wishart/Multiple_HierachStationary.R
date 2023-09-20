# Multiple Stationary Time Series simulated from perturbations with the same mean from a single time series
library(mvtnorm)
library(progress)
library(fda)
library(bayesplot)
library(ggplot2)
source("R/gradient_hierarch_Lambda.R")
source("R/posterior_hierarch_Lambda.R")
source("R/he_hierarch_Lambda.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")
source("R/Data_Generation/data_generation.R")
# Create time series 
# set parameters for generating data
# length of a single time series
n = 1000
# highest little j index value for the frequencies
J = floor((n-1) / 2)
# Frequency (\omega_j): defined on [0, 2\pi)
omega = (2 * pi * (0:J)) / n
# burn-in period for ARsim
burn = 50


# Need to Create ~ R copies of the time series and store it in a matrix
# Each column of the matrix contains a time series
# R : the number of independent stationary time series (R) 
R = 5

# call function to generate Krafty
data_gen = generate_Krafty(n = n, R = R)
#extract pieces from function
matrix_timeseries = data_gen$matrix_timeseries
theta_true = data_gen$theta_true

# Define Periodogram
# Define y_n(\omega_j) for the posterior function below
perio = (abs(mvfft(matrix_timeseries)) ^ 2 / n)

# subset perio for unique values, J = ceil((n-1) / 2) 
perio = perio[(0:J) + 1, , drop=FALSE]

par(mfrow = c(1,1))
plot(perio[,1], type = "l")

#################
# MCMC parameters
#################

# number of basis functions/number of beta values
B = 10
# Define lambda for the half-t prior: \pi(\tau^2 | \lambda) and \pi(\lambda)
lambda = 1
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

# prior variance for beta_0
sigmasquare = 100

# Set number of iterations
iter = 1000

# prior for Lambda : Wishart
# degrees of freedom
deg = B+1
# initial V
V = diag(B+1)

# Calculate ybar(omega_j)
y_bar = rowMeans(perio)



#######################
# Initialize parameters
#######################
# set tau^2 value
tausquared = 50
# The new D matrix that houses the prior variance of \beta^* 
Sigma = c(sigmasquare, D * tausquared)

# Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
# ncol: number of parameters (beta^*, tau^2)
# dim : (iter) x (B + 2)
Theta = matrix(NA, nrow = iter, ncol = B + 2)

# Create matrix of the basis functions
# fix fourier frequencies
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1
dim(Psi)
# Check X is orthogonal basis
round(crossprod(Psi),5)
# not orthogonal because we are not evaluating the periodogram at the full n-1 values.
# Initialize beta using least squares solution
# Using J amount of data for periodogram, can initialize beta this way:
betavalues = solve(crossprod(Psi), crossprod(Psi, log(y_bar)))

# Initialize bb_beta at the mean for the prior of bb_beta
bb_beta = tcrossprod(betavalues, rep(1,R))

# Initialize Lambda^{-1}
Lambda_inv = deg * V

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

# Create array to hold Lambda values
Lambda_array = array(data = NA, dim = c(iter, nrow(Lambda_inv), ncol(Lambda_inv)))
# Initalize first array with Lambda value
Lambda_array[1,,] = Lambda_inv

#####################
# MCMC Algorithm
#####################

#Rprof()
pb = progress_bar$new(total = iter - 1)
t1 = Sys.time()
for (g in 2:iter) {
  pb$tick()
  #g = 2
  #########################
  # tau^2 and lambda update
  #########################
  # Update \lambda and \tau with Gibbs Sampler
  lambda = 1/rgamma(1, (nu+1)/2, nu/tausquared + etasq)
  tausquared = 1/rgamma(1, (B + 1 + nu)/2, sum(betavalues[-1]^2 / D) / 2 + nu/lambda)
  # Update Theta matrix with new tau squared value
  Theta[g,B+2] = tausquared
  # Update Sigma with new tau^2 value
  Sigma = c(sigmasquare, D * tausquared)
  
  #####################
  # beta update
  #####################
  # Update \beta with Gibbs Sampler
  Lambda_solve = solve(diag(Sigma))
  betavalues = c(rmvnorm(n = 1, mean = solve(Lambda_solve + R * Lambda_inv) %*% Lambda_inv %*% bb_beta %*% rep(1,R),
                       sigma = solve(Lambda_solve + R * Lambda_inv)))
  #save new betavalue
  Theta[g, -(B+2)] = betavalues
  
  #####################
  # Lambda^{-1} update
  #####################
  # Update \Lambda^{-1} with Gibbs Sampler
  # obtain new Lambda_inv with Wishart distribution
  Lambda_inv = rWishart(n = 1, df = deg + R, Sigma = solve(tcrossprod(bb_beta - tcrossprod(betavalues, rep(1,R))) + diag(B+1)))[,,1]
  # save Lambda_inv 
  Lambda_array[g,,] = Lambda_inv
  
  ######################
  # bb_beta update :MH
  ######################
  # Update \mathbb{B} with Metropolis Hastings Sampler
  for(r in 1:R){
    #r = 1
    br = bb_beta[,r]
    # Maximum A Posteriori (MAP) estimate : finds the \beta that gives us the mode of the conditional posterior of \beta conditioned on y
    map <- optim(par = br, fn = posterior_hierarch_Lambda, gr = gradient_hierarch_Lambda, method ="BFGS", control = list(fnscale = -1),
                 Psi = Psi, sumPsi = sumPsi, y = perio[,r], b = betavalues, lambda = Lambda_inv)$par
    # Call the hessian function
    norm_precision <- he_hierarch_Lambda(br = map, Psi = Psi, y = perio[,r], lambda = Lambda_inv) * -1
    # Calculate the \beta^* proposal, using Cholesky Sampling
    betaprop <- Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
    # Calculate acceptance ratio
    prop_ratio <- min(1, exp(posterior_hierarch_Lambda(br = betaprop, b = betavalues, Psi = Psi, sumPsi = sumPsi, y = perio[,r], lambda = Lambda_inv) -
                               posterior_hierarch_Lambda(br = br, b = betavalues, Psi = Psi, sumPsi = sumPsi, y = perio[,r], lambda = Lambda_inv)))
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
# Rprof(NULL)
# summaryRprof()


#######################
# Plots and Diagnostics
#######################
# Need to update this code for the hierarchical set-up
# Remove burn-in
burnin = 10
Theta = Theta[-(1:burnin),]
bb_beta_array = bb_beta_array[-(1:burnin),,]
Lambda_array = Lambda_array[-(1:burnin),,]
r = 5


# Create spectral densities for the first series
specdens = exp(Psi %*% t(bb_beta_array[,,r]))
par(mfrow = c(1, 1))
plot(x =c(), y=c(), xlim = c(0,3), ylim = range(log(specdens)), ylab = "Spectral Density", xlab = "omega",
     main = "Spectral Density Estimates \nwith True Spectral Density")
for(h in sample(ncol(specdens), 100, replace = FALSE)){
  lines(x = omega, y = log(specdens[,h]), col = rgb(0, 0, 0, 0.2))
}
lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_vec[r])), col = "red", lwd = 2)
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Estimate", "True"))


theta_vec

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



















