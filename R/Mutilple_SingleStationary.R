# Multiple Single Stationary Time Series from the same underlying spectra
library(mvtnorm)
library(progress)
set.seed(22)
source("R/posterior_multiple.R")
source("R/gr_multiple.R")
source("R/he_multiple.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")

# Create a single time series
# set hyper-parameters
# length of a single time series
n <-  1000
# Frequency (\omega_j): defined on [0, 2\pi)
# for j = 0,...,n-1
omega <- (2 * pi * (0:(n-1)))/n
# burn-in period for ARsim
burn <- 50

# Pick an AR() process: AR(1), AR(2), AR(3)
# Create coefficient \phi
# For AR(1)
phi <- 0.5

# For AR(2)
# phi <- c(1.4256, -0.9)

# For AR(3)
# phi <- c(1.4256, -0.7344, 0.1296)

# Need to Create ~ B copies of the time series and store it in a matrix
# Each column of the matrix contains a time series
# B : the number of independent stationary time series (B) 
B <- 500
# create matrix to store the time series: (B x n)
matrix_timeseries <- matrix(NA, nrow = n, ncol = B)
for(r in 1:B){
  matrix_timeseries[,r] <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)
}


# Plot the time series that will be used
par(mfrow = c(5,2))
for(i in 1:(ncol(matrix_timeseries))){
  plot(x = matrix_timeseries[,i], type = "l")
}

# Define Periodogram
# Define y_n(\omega_j) for the posterior function below
perio <- (abs(mvfft(matrix_timeseries)) ^ 2 / n)

par(mfrow = c(5,2))
for(i in 1:(ncol(perio))){
  plot(omega, perio[,i], type = "l", xlim = c(0,pi))
}

# Calculate ybar(omega_j)
y_bar <- rowMeans(perio)

#################
# MCMC parameters
#################

# number of basis functions/number of beta values
K <- 10
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
# D = rep(1, K)

# Yakun's D
D = c((sqrt(2)*pi*(1:K))^(-2))

# exponential decay D
# D = exp(0.12 * -(0:(K-1)))

# prior variance for beta_0
sigmasquare = 100

# Set number of iterations
iter = 10000

#######################
# Initialize parameters
#######################
# set tau^2 value
tausquared <- 50
# The new D matrix that houses the prior variance of \beta^* 
Sigma = c(sigmasquare, D * tausquared)

# Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
# ncol: number of parameters (beta^*, tau^2)
# dim : (iter) x (K + 2)
Theta <- matrix(NA, nrow = iter, ncol = K + 2)

# Create matrix of the basis functions
# fix fourier frequencies
X <- outer(X = omega, Y = 0:K, FUN = function(x,y){sqrt(2) * cos(y * x)})
# redefine the first column to be 1's
X[,1] <- 1
# Initialize beta using least squares solution
betavalues <- c(crossprod(X,log(y_bar))) / n

# Check X is orthogonal basis
round(crossprod(X),5)
# Specify Sum of X for the posterior function later
# 1^T_n X part in the paper: identical to colSums but is a faster calculation
sumX <- c(crossprod(rep(1, nrow(X)), X))
# Initialize first row of Theta
Theta[1,] <- c(betavalues, tausquared)

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
  b = Theta[g - 1, 1:(K+1)]
  # tau^squared of most recent iteration:
  tsq = Theta[g - 1, K + 2]
  ##########################
  # Metropolis Hastings Step
  ##########################
  # Maximum A Posteriori (MAP) estimate : finds the \beta^* that gives us the mode of the conditional posterior of \beta^* conditioned on y
  map <- optim(par = b, fn = posterior_multiple, gr = gr_multiple, method ="BFGS", control = list(fnscale = -1),
               X = X, sumX = sumX, y_bar = y_bar, D = Sigma, B = B)$par
  # Call the hessian function
  norm_precision <- he_multiple(b = map, X = X, sumX = sumX, y_bar = y_bar, D = Sigma , B = B)
  # Calculate the \beta^* proposal, using Cholesky Sampling
  betaprop <- Chol_sampling(Lt = chol(norm_precision), d = K + 1, beta_c = map)
  # Calculate acceptance ratio
  prop_ratio <- min(1, exp(posterior_multiple(b = betaprop, X = X, sumX = sumX, y_bar = y_bar,  D = Sigma, B = B) -
                             posterior_multiple(b = b, X = X, sumX = sumX, y_bar = y_bar,  D = Sigma, B = B)))
  # Create acceptance decision
  accept <- runif(1)
  if(accept < prop_ratio){
    # Accept betaprop as new beta^*
    Theta[g, -(K+2)] <- betaprop
  }else{
    # Reject betaprop as new beta^*
    Theta[g, -(K+2)] <- b
  }
  ##############################
  # Tau^2 Update: Gibbs Sampler: conditional conjugate prior for the half-t
  ##############################
  # 1/rgamma is the same as invgamma (so we don't need the other library)
  lambda = 1/rgamma(1, (nu+1)/2, nu/tsq + etasq)
  newtsq = 1/rgamma(1, (K + 1 + nu)/2, sum(Theta[g, -(K+2)]^2 / D) / 2 + nu/lambda)
  # Update Theta matrix with new tau squared value
  Theta[g,K+2] <- newtsq
  # Update Sigma with new tau^2 value
  Sigma = c(sigmasquare, D * newtsq)
}
#Rprof(NULL)
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
specdens <- exp(X %*% t(Theta[ ,-(K+2)]))
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
pdf(file = "Posterior_Mean_Mulitple.pdf",
    width = 10,
    height = 5,)
par(mfrow = c(1, 1))
plot(x =c(), y=c(), xlim = c(0,pi), ylim = range(log(specdens)), ylab = "Spectral Density", xlab = "omega",
     main = "Posterior Mean and\n 95% Credible Interval")
polygon(x = c(omega,rev(omega)), y = log(c(summary_stats$lower, rev(summary_stats$upper))), col = "darkgrey", border = NA)
#lines(x = omega, y = summary_stats$lower, lty = 2, col = "darkgrey")
lines(x = omega, y = log(summary_stats$mean), col = "black")
#lines(x = omega, y = summary_stats$upper, lty = 2, col = "darkgrey")
lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "red", lwd = 2)
#lines(x = omega, y = exp(perio), col = "lightgrey")
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Posterior Mean", "True Spectral Density"))
#lines(y = log(y_bar), x = omega, col = "orange")
dev.off()

mean((arma_spec(omega = omega, phi = phi) - summary_stats$mean)^2)
# for n = 1000
# 0.0211


##############
# Trace Plots
##############

# Beta Trace Plots
pdf(file = "BetaTrace_AR1_Multiple.pdf",
    width = 12,
    height = 6)
par(mfrow = c(4,4), mar = c(4.2, 4.2, 2, 0.2))
for(m in 1:(K + 1)){
  plot(Theta[,m], type = "l")
}
mtext(text = "Beta Trace Plots AR(1)", outer = TRUE, line = -1.5)
dev.off()

# plot tau^2 estimate
pdf(file = "tausquaredTrace_AR1_Multiple.pdf",
    width = 12,
    height = 6)
par(mfrow = c(1,1))
plot(Theta[,K+2], type = "l")
dev.off()

mean(abs(sign(diff(Theta[,1]))))
# 25% for S = 500; n = 1000; K = 10 

# Check how many point-wise times the red line is outside the 95th credible interval
mean(arma_spec(omega = omega, phi = phi) > summary_stats$upper | arma_spec(omega = omega, phi = phi) < summary_stats$lower)
