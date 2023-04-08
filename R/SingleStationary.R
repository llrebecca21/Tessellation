# Single Stationary Time Series

set.seed(1080)

# Create single time series: AR(1)

# set hyper-parameters

# length of time series
n <- 1000
# burn-in period
burn <- 50
# Create coefficient phi
phi <- 0.5

# Model 1 : AR(1) with phi = 0.5
ts1 <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)

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
omega <- (0:(n-1)) / n

# Set number of iterations
iter <- 100
  
#######################
# Initialize parameters
#######################
# set tau^2 value
tausquared <- 10
# set intercept term
alpha0 <- 0
# set beta
betavalues <- rep(0, K)

# Create matrix to store samples
# ncol: number of parameters (alpha0, Beta, tau^2)
Theta <- matrix(NA, nrow = iter, ncol = K + 2)

# Create matrix of the basis functions
# unscale by n (i.e. 2/n)
# fix fourier frequencies by multiplying by 2 (inside cosine function)
X <- outer(X = omega, Y = 1:K, FUN = function(x,y){sqrt(2/n)*cos(2 * pi * y * x)})
# Check X is orthonormal
round(crossprod(X),5)

# Initialize first row of Theta
Theta[1,] <- c(alpha0, betavalues, tausquared)

#####################
# MCMC Algorithm
#####################

for (i in 2:iter) {
  # Metropolis Hastings Step
  
  
  # Tau Update: Gibbs Sampler: Inverse CDF sampler
  # truncated gamma
  # draw u first: corresponds to a valid Gamma CDF value
  u <- runif(1, min = pgamma(q = 1/maxtausquared, shape = K/2 - 1, rate = 1/2 * crossprod(Theta[i,-c(1, K+2)])) , max = 1)
  # Recovering corresponding inverse tau squared
  invtausquarednew <- qgamma(p = u, shape = K/2 -1, rate = 1/2 * crossprod(Theta[i, -c(1,K+2)]))
  # invert to get new tau squared value
  newtau <- 1/invtausquarednew
  # Update Theta matrix with new tau squared value
  Theta[i,K+2] <- newtau
}





