# Multiple Single Stationary Time Series from the same underlying spectra
library(mvtnorm)
library(progress)
set.seed(1080)
source("R/posterior_multiple.R")
source("R/gr_multiple.R")
source("R/he_multiple.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")

# Create single time series: AR(1)

# set hyper-parameters
# length of time series
n <-  1000
# burn-in period for ARsim
burn <- 50

# Create coefficient phi
# For AR(1)
phi <- 0.5

# AR(2)
# phi <- c(1.4256, -0.9)

# For AR(3)
# phi <- c(1.4256, -0.7344, 0.1296)

# Need to Create ~ 10 copies of the time series and store it in a matrix
# Each column of the matrix contains a time series
num_timeseries <- 10
# create matrix to store the time series
matrix_timeseries <- matrix(NA, nrow = n, ncol = num_timeseries)
for(r in 1:num_timeseries){
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
K <- 11
# alphabeta stores number of beta values (beta_{1:K}) + intercept term (alpha_0)
# alphabeta <- K + 1
# prior intercept variance (variance associated with the alpha_0 prior)
# sigmasalpha <- 100
# Define omega for the Basis Functions
omega <- (2 * pi * (0:(n-1)))/n
# Define lambda for the half-t prior
lambda = 1
# Define degrees of freedom for half-t prior
# Cauchy nu = 1
nu = 3
# Define eta as the other scale parameter for half-t prior
# eta = 1 gives standard Cauchy; higher eta gives wider Cauchy
eta = 1
# Define D's main diagonal
# D is a measure of prior variance
# Identity D:
# D = rep(1, K)

# D = c((sqrt(2)*pi*(1:K))^(-2)) * 100
D = c(1, (sqrt(2)*pi*(1:(K-1)))^(-2))

# D = exp(0.12 * -(0:(K-1)))

# Set number of iterations
iter <- 10000

#######################
# Initialize parameters
#######################
# set tau^2 value
tausquared <- 50
# set intercept term
# alpha0 <- 0
# set beta
# betavalues <- rnorm(K,mean = 0, sd = sqrt(tausquared))
# Initialize beta at 0
# betavalues <- rep(0,K)

# Initialize beta using least squares solution
betavalues <- c(crossprod(X,log(y_bar))) / c(n, rep(n/2, K-1))

# Create matrix to store samples
# ncol: number of parameters (Beta, tau^2)
Theta <- matrix(NA, nrow = iter, ncol = K + 1)

# Create matrix of the basis functions
# unscale by n (i.e. 2/n)
# fix fourier frequencies
X <- outer(X = omega, Y = 0:(K-1), FUN = function(x,y){cos(y * x)})
# Check X is orthogonal basis
dim(X)
# Fix the scaling of the first column
#X[,1] <- X[,1]/sqrt(2)
round(crossprod(X),5)
# Specify Sum of X for the posterior function later
# 1^T_n X Beta part in the paper (excluding the Beta)
sumX <- c(crossprod(rep(1, nrow(X)), X))
length(sumX)
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
  # Extract beta and tau^2 from theta
  # beta:
  b = Theta[g - 1, 1:K]
  # tau^squared:
  tsq = Theta[g - 1, K + 1]
  # Metropolis Hastings Step
  # Maximum A Posteriori (MAP) estimate : finds the alpha and beta that gives us the mode of the conditional posterior of beta and alpha_0 conditioned on y
  map <- optim(par = b, fn = posterior_multiple, gr = gr_multiple, method ="BFGS", control = list(fnscale = -1),
               X = X, sumX = sumX, tsq = tsq, y_bar = y_bar, D = D, num_timeseries = num_timeseries)$par
  # Call the hessian function
  norm_precision <- he_multiple(b = map, X = X, sumX = sumX, tsq = tsq, y_bar = y_bar, D = D , num_timeseries = num_timeseries)
  # Calculate the alpha beta proposals, using Cholesky Sampling
  betaprop <- Chol_sampling(Lt = chol(norm_precision), d = K, beta_c = map)

  # calculate acceptance ratio
  prop_ratio <- min(1, exp(posterior_multiple(b = betaprop, X = X, sumX = sumX, tsq = tsq, y_bar = y_bar,  D = D, num_timeseries = num_timeseries) -
                             posterior_multiple(b = b, X = X, sumX = sumX, tsq = tsq, y_bar = y_bar,  D = D, num_timeseries = num_timeseries)))
  # create acceptance decision
  accept <- runif(1)
  if(accept < prop_ratio){
    # accept betaprop as new beta
    Theta[g, -(K+1)] <- betaprop
  }else{
    Theta[g, -(K+1)] <- b
  }
  # Tau^2 Update: Gibbs Sampler: conditional conjugate prior for the half-t
  # 1/gamma is the same as invgamma (so we dont need the other library)
  lambda = 1/rgamma(1, (nu+1)/2, nu/tsq + eta^2)
  newtsq = 1/rgamma(1, (K + nu)/2, sum(Theta[g, -(K+1)]^2 / D) / 2 + nu/lambda)
  # Update Theta matrix with new tau squared value
  Theta[g,K+1] <- newtsq
  #Theta[g,K+2] <- 100
}
#Rprof(NULL)
# summaryRprof()

# Remove burn-in
burnin <- 10
Theta <- Theta[-(1:burnin),]

# Plot the Spectral Density Estimates
#pdf(file = "Spectral_Density_Estimates.pdf",
#    width = 10,
#    height = 5,)
specdens <- exp(X %*% t(Theta[ ,-(K+1)]))
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
plot(x =c(), y=c(), xlim = c(0,pi), ylim = range(specdens), ylab = "Spectral Density", xlab = "omega",
     main = "Posterior Mean and\n 95% Credible Interval")
polygon(x = c(omega,rev(omega)), y = c(summary_stats$lower, rev(summary_stats$upper)), col = "darkgrey", border = NA)
#lines(x = omega, y = summary_stats$lower, lty = 2, col = "darkgrey")
lines(x = omega, y = summary_stats$mean, col = "black")
#lines(x = omega, y = summary_stats$upper, lty = 2, col = "darkgrey")
lines(x = omega, y = arma_spec(omega = omega, phi = phi), col = "red", lwd = 2)
#lines(x = omega, y = exp(perio), col = "lightgrey")
legend("topright", col = c("black", "red"), lwd = c(1,2), legend = c("Posterior Mean", "True Spectral Density"))
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
par(mfrow = c(3,4), mar = c(4.2, 4.2, 2, 0.2))
for(m in 1:(K)){
  plot(Theta[,m], type = "l")
}
mtext(text = "Beta Trace Plots AR(1)", outer = TRUE, line = -1.5)
dev.off()

# plot tau^2 estimate
pdf(file = "tausquaredTrace_AR1_Multiple.pdf",
    width = 12,
    height = 6)
par(mfrow = c(1,1))
plot(Theta[,K+1], type = "l")
dev.off()

mean(abs(sign(diff(Theta[,1]))))
# 27% acceptance rate approximately



