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

# 31 : perfect at edges use for example that this is a seed thing


# Create a single time series
# set hyper-parameters
# length of a single time series
n = 1000
# highest little j index value for the frequencies
J = floor((n-1) / 2)
# Frequency (\omega_j): defined on [0, 2\pi)
# for j = 0,...,n-1
# omega = (2 * pi * (0:(n-1)))/n
omega = (2 * pi * (0:J)) / n
# burn-in period for ARsim
burn = 50

# Pick an AR() process: AR(1), AR(2), AR(3)
# Create coefficient \phi
# For AR(1)
phi = 0.5

# For AR(2)
#phi <- c(1.4256, -0.9)

# For AR(3)
#phi <- c(1.4256, -0.7344, 0.1296)

# Need to Create ~ R copies of the time series and store it in a matrix
# Each column of the matrix contains a time series
# R : the number of independent stationary time series (R) 
R = 5
# create matrix to store the time series: (R x n)
matrix_timeseries = matrix(NA, nrow = n, ncol = R)
for(r in 1:R){
  matrix_timeseries[,r] <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)
}
dim(matrix_timeseries)
plot(matrix_timeseries[,1], type = "l")

# Plot the time series that will be used
# par(mfrow = c(5,2))
# for(i in 1:(ncol(matrix_timeseries))){
#   plot(x = matrix_timeseries[,i], type = "l")
# }

# Define Periodogram
# Define y_n(\omega_j) for the posterior function below
perio = (abs(mvfft(matrix_timeseries)) ^ 2 / n)
dim(perio)
# subset perio for unique values, J = ceil((n-1) / 2) 
perio = perio[(0:J) + 1, , drop=FALSE]
dim(perio)
par(mfrow = c(1,1))
plot(perio[,1], type = "l")

# par(mfrow = c(5,2))
# for(i in 1:(ncol(perio))){
#   plot(omega, perio[,i], type = "l", xlim = c(0,pi))
# }

# Calculate ybar(omega_j)
y_bar = rowMeans(perio)

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

# Define D's main diagonal : Choose either identity, Yakun's D, or exponential decay D
# D is a measure of prior variance for \beta_1 through \beta_K

# Identity D:
# D = rep(1, B)

# Rebecca's D
D = 1 / (4 * pi * (1:B)^2)

# exponential decay D
# D = exp(0.12 * -(0:(B-1)))

# prior variance for beta_0
sigmasquare = 100

# Set number of iterations
iter = 1000

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
# Using the full data n-1 for periodogram can use this to initialize beta:
#betavalues = c(crossprod(Psi,log(y_bar))) / n

# Using J amount of data for periodogram, can initialize beta this way:
betavalues = solve(crossprod(Psi), crossprod(Psi, log(y_bar)))


# Specify Sum of X for the posterior function later
# Specify Sum of X for the posterior function later
# 1^T_n X part in the paper: identical to colSums but is a faster calculation
sumPsi = c(crossprod(rep(1, nrow(Psi)), Psi))
# Initialize first row of Theta
Theta[1,] = c(betavalues, tausquared)

#####################
# MCMC Algorithm
#####################

#Rprof()
pb = progress_bar$new(total = iter - 1)
for (g in 2:iter) {
  pb$tick()
  # g = 2
  # Extract \beta^* and tau^2 from theta
  # beta^* of most recent iteration:
  b = Theta[g - 1, 1:(B+1)]
  # tau^squared of most recent iteration:
  tsq = Theta[g - 1, B + 2]
  ##########################
  # Metropolis Hastings Step
  ##########################
  # Maximum A Posteriori (MAP) estimate : finds the \beta^* that gives us the mode of the conditional posterior of \beta^* conditioned on y
  map <- optim(par = b, fn = posterior_multiple, gr = gr_multiple, method ="BFGS", control = list(fnscale = -1),
               Psi = Psi, sumPsi = sumPsi, y_bar = y_bar, D = Sigma, R = R)$par
  # Call the hessian function
  norm_precision <- he_multiple(b = map, Psi = Psi, y_bar = y_bar, D = Sigma , R = R)
  # Calculate the \beta^* proposal, using Cholesky Sampling
  betaprop <- Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
  # Calculate acceptance ratio
  prop_ratio <- min(1, exp(posterior_multiple(b = betaprop, Psi = Psi, sumPsi = sumPsi, y_bar = y_bar,  D = Sigma, R = R) -
                             posterior_multiple(b = b, Psi = Psi, sumPsi = sumPsi, y_bar = y_bar,  D = Sigma, R = R)))
  # Create acceptance decision
  accept <- runif(1)
  if(accept < prop_ratio){
    # Accept betaprop as new beta^*
    Theta[g, -(B+2)] <- betaprop
  }else{
    # Reject betaprop as new beta^*
    Theta[g, -(B+2)] <- b
  }
  ##############################
  # Tau^2 Update: Gibbs Sampler: conditional conjugate prior for the half-t
  ##############################
  # 1/rgamma is the same as invgamma (so we don't need the other library)
  lambda = 1/rgamma(1, (nu+1)/2, nu/tsq + etasq)
  newtsq = 1/rgamma(1, (B + 1 + nu)/2, sum(Theta[g, -c(1,B+2)]^2 / D) / 2 + nu/lambda)
  # Update Theta matrix with new tau squared value
  Theta[g,B+2] = newtsq
  # Update Sigma with new tau^2 value
  Sigma = c(sigmasquare, D * newtsq)
}
# Rprof(NULL)
# summaryRprof()


#######################
# Plots and Diagnostics
#######################
# Remove burn-in
burnin <- 10
Theta <- Theta[-(1:burnin),]

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


# Plot with the bounds:
# pdf(file = "Posterior_Mean_Mulitple.pdf",
#     width = 10,
#     height = 5,)
par(mfrow = c(1, 1), mar = c(4,4,4,1) + .1)
plot(x =c(), y=c(), xlim = c(0,pi), ylim = range(log(specdens)), ylab = "Spectral Density", xlab = "omega",
     main = sprintf("Posterior Mean and\n 95%s Credible Interval, B = %g, R = %g, n = %g", "%", B, R, n))
polygon(x = c(omega,rev(omega)), y = log(c(summary_stats$lower, rev(summary_stats$upper))), col = "darkgrey", border = NA)
#lines(x = omega, y = summary_stats$lower, lty = 2, col = "darkgrey")
lines(x = omega, y = log(summary_stats$mean), col = "black")
#lines(x = omega, y = summary_stats$upper, lty = 2, col = "darkgrey")
lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "red", lwd = 2)
#lines(x = omega, y = exp(perio), col = "lightgrey")
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Posterior Mean", "True Spectral Density"))
#lines(y = log(y_bar), x = omega, col = "orange")
#dev.off()

mean((arma_spec(omega = omega, phi = phi) - summary_stats$mean)^2)
# for n = 1000
# 0.0211


##############
# Trace Plots
##############

# Beta Trace Plots
#pdf(file = "BetaTrace_AR1_Multiple.pdf",
#    width = 12,
#    height = 6)
par(mfrow = c(4,4), mar = c(4.2, 4.2, 2, 0.2))
for(m in 1:(B + 1)){
  plot(Theta[,m], type = "l")
}
mtext(text = "Beta Trace Plots AR(1)", outer = TRUE, line = -1.5)
#dev.off()

# plot tau^2 estimate
#pdf(file = "tausquaredTrace_AR1_Multiple.pdf",
#    width = 12,
#    height = 6)
par(mfrow = c(1,1))
plot(Theta[,B+2], type = "l")
#dev.off()

mean(abs(sign(diff(Theta[,1]))))
# seed = 22
# 25% for B = 500 ; n = 1000; K = 10 
# 26% for B = 1000; n = 1000; K = 10
# 27% for B = 1   ; n = 1000; K = 10
# 42% for B = 1   ; n = 1000; K = 5

# Check how many point-wise times the red line (True Spectral Density) is outside the 95th credible interval
mean(arma_spec(omega = omega, phi = phi) > summary_stats$upper | arma_spec(omega = omega, phi = phi) < summary_stats$lower)
# seed = 22
# 0.038 for B = 500 ; n = 1000; K = 10 
# 0.162 for B = 1000; n = 1000; K = 10
# 0.161 for B = 1   ; n = 1000; K = 10
# 0.187 for B = 1   ; n = 1000; K = 5




