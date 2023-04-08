# Single Stationary Time Series
library(mvtnorm)
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
omega <- (1:ceiling(n/2) - 1) / n

# Set number of iterations
iter <- 10000
  
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
X <- outer(X = omega, Y = 1:K, FUN = function(x,y){sqrt(2/length(omega))*cos(4 * pi * y * x)})
# Check X is orthonormal basis
round(crossprod(X),5)
# Specify Sum of X for the posterior function later
# 1^T_n X Beta part in the paper (excluding the Beta)
sumX <- crossprod(rep(1, nrow(X)), X)

# Initialize first row of Theta
Theta[1,] <- c(alpha0, betavalues, tausquared)

#####################
# MCMC Algorithm
#####################
# Define y_n(\omega_j) for the posterior function below
perio <- log((abs(fft(ts1)) ^ 2 / n)[1:ceiling(n/2)])

plot(perio, type = "l")
length((abs(fft(ts1)) ^ 2 / n))

# Define Posterior Function for Beta and alpha_0
post_func <- function(b, a, t){
  exp(-1/2 * (crossprod(b) / t + a^2 / sigmasalpha + nrow(X) * a + sumX %*% b +
                sum(exp(perio - a - X %*% b))))
}


for (i in 2:iter) {
  # Metropolis Hastings Step
  betaprop <- rmvnorm(n = 1, mean = Theta[i - 1, -(K+2)], sigma = 0.03 * diag(K+1))
  # calculate acceptance ratio
  prop_ratio <- min(1, post_func(b = betaprop[-1], a = betaprop[1], t = Theta[i - 1, K + 2]) / post_func(b = Theta[i-1, -c(1, K+2)], a = Theta[i-1, 1], t = Theta[i - 1, K + 2]) )
  # create acceptance decision
  accept <- rbinom(1, 1, prop_ratio)
  if(accept == 1){
    # accept betaprop as new alpha and beta
    Theta[i, -(K+2)] <- betaprop
  }else{
    Theta[i, -(K+2)] <- Theta[i - 1, -(K+2)]
  }
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

Theta <- Theta[-(1:500),]

plot(Theta[,1], type = "l")
plot(Theta[,2], type = "l")
plot(Theta[,K+2], type = "l")

