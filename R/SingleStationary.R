# Single Stationary Time Series
library(mvtnorm)
set.seed(1080)
source("R/whittle_post.R")
source("R/gr_single.R")
source("R/he_single.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")

# Create single time series: AR(1)

# set hyper-parameters
# length of time series
n <-  2000
# burn-in period
burn <- 50
# Create coefficient phi
phi <- 0.5

# Model 1 : AR(1) with phi = 0.5
# Rprof()
ts1 <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)
# Rprof(NULL)
# summaryRprof()
#################
# MCMC parameters
#################

# number of basis functions/number of beta values
K <- 10
# nbeta stores number of beta values (beta_{1:K}) + intercept term (alpha_0)
alphabeta <- K + 1
# prior intercept variance (variance associated with the alpha_0 prior)
sigmasalpha <- 100
# maximum value for tau^2 (Indicator in paper).
maxtausquared <- 1000
# Define omega for the Basis Functions
omega <- (2 * pi * (0:(n-1)))/n
# Define lambda for the half-t prior
lambda = 1
# Define degrees of freedom for half-t prior
# Cauchy
nu = 1 
# Define eta as the other scale parameter for half-t prior
# eta = 1 gives standard Cauchy; higher eta gives wider Cauchy
eta = 1
# Define D's main diagonal
# Identity D:
# D = rep(1, K)
# Slow Decaying D
D = exp(0.12 * -(0:(K-1)))
plot(D)
# Set number of iterations
iter <- 10000
  
#######################
# Initialize parameters
#######################
# set tau^2 value
tausquared <- 0.1
# set intercept term
alpha0 <- 0
# set beta
betavalues <- rnorm(K,mean = 0, sd = sqrt(tausquared))

# Create matrix to store samples
# ncol: number of parameters (alpha0, Beta, tau^2)
Theta <- matrix(NA, nrow = iter, ncol = K + 2)

# Create matrix of the basis functions
# unscale by n (i.e. 2/n)
# fix fourier frequencies by multiplying by 2 (inside cosine function)
X <- outer(X = omega, Y = 1:K, FUN = function(x,y){cos(y * x) / sqrt(n / 2)})
# Check X is orthonormal basis
dim(X)
# 2000 rows : n
# 20: K
round(crossprod(X),5)
# Specify Sum of X for the posterior function later
# 1^T_n X Beta part in the paper (excluding the Beta)
sumX <- c(crossprod(rep(1, nrow(X)), X))
length(sumX)
# Initialize first row of Theta
Theta[1,] <- c(alpha0, betavalues, tausquared)

#####################
# MCMC Algorithm
#####################
# Define y_n(\omega_j) for the posterior function below

perio <- log((abs(fft(ts1)) ^ 2 / n))

plot(omega, perio, type = "l")
length((abs(fft(ts1)) ^ 2 / n))

Rprof()
for (g in 2:iter) {
  # g = 2
  # Extract alpha, beta and tau^2 from theta
  # alpha:
  a = Theta[g - 1, 1]
  # beta:
  b = Theta[g - 1, -c(1, K+2)]
  # tau^squared:
  tsq = Theta[g - 1, K + 2]
  # Metropolis Hastings Step
  # Maximum A Posteriori (MAP) estimate : finds the alpha and beta that gives us the mode of the conditional posterior of beta and alpha_0 conditioned on y
  map <- optim(par = c(a,b), fn = whittle_post, gr = gr_single, method ="BFGS", control = list(fnscale = -1),
               X = X, sumX = sumX, tsq = tsq, perio = perio, sigmasalpha = sigmasalpha, D = D)$par
  # Call the hessian function
  norm_precision <- he_single(ab = map, X = X, sumX = sumX, tsq = tsq, perio = perio, sigmasalpha = sigmasalpha, D = D)
  # Calculate the alpha beta proposals, using Cholesky Sampling
  betaprop <- Chol_sampling(Lt = chol(norm_precision), d = K + 1, beta_c = map)
  
  #betaprop <- rmvnorm(n = 1, mean = Theta[g - 1, -(K+2)], sigma = 0.03 * diag(K+1))
  # calculate acceptance ratio
  prop_ratio <- min(1, exp(whittle_post(ab = betaprop, X = X, sumX = sumX, tsq = tsq, perio = perio, sigmasalpha = sigmasalpha, D = D) -
                      whittle_post(ab = c(a,b), X = X, sumX = sumX, tsq = tsq, perio = perio, sigmasalpha = sigmasalpha, D = D)))
  # create acceptance decision
  accept <- runif(1)
  if(accept < prop_ratio){
    # accept betaprop as new alpha and beta
    Theta[g, -(K+2)] <- betaprop
  }else{
    Theta[g, -(K+2)] <- c(a,b)
  }
  # Tau^squared Update: Gibbs Sampler: conditional conjugate prior for the half-t
  # 1/gamma is the same as invgamma (so we dont need the other library)
  lambda = 1/rgamma(1, (nu+1)/2, nu/tsq + eta^2)
  newtsq = 1/rgamma(1, (K + nu)/2, crossprod(Theta[g,-c(1,K+2)]) / 2 + nu/lambda)
  # Update Theta matrix with new tau squared value
  Theta[g,K+2] <- newtsq
}
Rprof(NULL)
summaryRprof()

# Remove burn-in
Theta <- Theta[-(1:burn),]

plot(Theta[,1], type = "l")
# betas
plot(Theta[,2], type = "l")


par(mfrow = c(2,5), mar = c(4.2, 4.2, 2, 0.2))
for(m in 2:(K+1)){
  plot(Theta[,m], type = "l")
}

# plot tau estimate
plot(Theta[,K+2], type = "l")

# Plot the Spectral Density Estimates
pdf(file = "Spectral_Density_Estimates.pdf",
    width = 10,
    height = 5,)
specdens <- exp(cbind(1,X) %*% t(Theta[ ,-(K+2)]))
par(mfrow = c(1, 1))
plot(x =c(), y=c(), xlim = range(omega), ylim = range(specdens), ylab = "Spectral Density", xlab = "omega",
     main = "Spectral Density Estimates \nwith True Spectral Density")
for(h in sample(ncol(specdens), 1000, replace = FALSE)){
  lines(x = omega, y = specdens[,h], col = rgb(0, 0, 0, 0.2))
}
lines(x = omega, y = arma_spec(omega = omega, phi = phi), col = "red", lwd = 2)
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Estimate", "True"))
dev.off()

dim(specdens)
# n X (numiter - burn)

# Create Data frame to store the lower bound and upper bound and mean
summary_stats <- data.frame("lower" = apply(specdens, 1, FUN = function(x){quantile(x, .025)}), "mean" = rowMeans(specdens),
                            "upper" = apply(specdens, 1, FUN = function(x){quantile(x, 0.975)}))


# Plot with the bounds:
pdf(file = "Posterior_Mean.pdf",
    width = 10,
    height = 5,)
par(mfrow = c(1, 1))
plot(x =c(), y=c(), xlim = range(omega), ylim = range(specdens), ylab = "Spectral Density", xlab = "omega",
     main = "Posterior Mean and\n 95% Confidence Interval")
polygon(x = c(omega,rev(omega)), y = c(summary_stats$lower, rev(summary_stats$upper)), col = "darkgrey", border = NA)
#lines(x = omega, y = summary_stats$lower, lty = 2, col = "darkgrey")
lines(x = omega, y = summary_stats$mean, col = "black")
#lines(x = omega, y = summary_stats$upper, lty = 2, col = "darkgrey")
lines(x = omega, y = arma_spec(omega = omega, phi = phi), col = "red", lwd = 2)
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Posterior Mean", "True Spectral Density"))
dev.off()

# Metropolis Hastings Step
# Bring in the Gaussian Approximation with BFGS Optimization
# use log to write the acceptance step
# change index i to g

# Tau Update
# Change update to half-t distribution

##### Notes: 

# get the derivation for the f(\omega) ~ N( , \tau^2 D) (choice of D)
# Whittle Posterior : Do not use the Identity, use the D
# Whittle, Gradient and the Hessian, need the prior













