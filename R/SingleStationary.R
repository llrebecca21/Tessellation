# Single Stationary Time Series
library(mvtnorm)
set.seed(1080)
source("R/whittle_post.R")
source("R/gr_single.R")

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
sumX <- c(crossprod(rep(1, nrow(X)), X))
length(sumX)
# Initialize first row of Theta
Theta[1,] <- c(alpha0, betavalues, tausquared)

#####################
# MCMC Algorithm
#####################
# Define y_n(\omega_j) for the posterior function below
perio <- log((abs(fft(ts1)) ^ 2 / n)[1:ceiling(n/2)])

plot(perio, type = "l")
length((abs(fft(ts1)) ^ 2 / n))

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
  # Maximum A Posteriori (MAP) estimate : finds the alpha and beta that gives us the mode of the conditional posterior
  map <- optim(par = c(a,b), fn = whittle_post, gr = gr_single, method ="BFGS", control = list(fnscale = -1),
               X = X, sumX = sumX, tsq = tsq, perio = perio, sigmasalpha = sigmasalpha)
  
  # 
  
  #betaprop <- rmvnorm(n = 1, mean = Theta[g - 1, -(K+2)], sigma = 0.03 * diag(K+1))
  # calculate acceptance ratio
  prop_ratio <- min(1, post_func(b = betaprop[-1], a = betaprop[1], t = Theta[i - 1, K + 2]) / post_func(b = Theta[g-1, -c(1, K+2)], a = Theta[g-1, 1], t = Theta[g - 1, K + 2]) )
  # create acceptance decision
  accept <- rbinom(1, 1, prop_ratio)
  if(accept == 1){
    # accept betaprop as new alpha and beta
    Theta[g, -(K+2)] <- betaprop
  }else{
    Theta[g, -(K+2)] <- Theta[g - 1, -(K+2)]
  }
  # Tau Update: Gibbs Sampler: Inverse CDF sampler
  # truncated gamma
  # draw u first: corresponds to a valid Gamma CDF value
  u <- runif(1, min = pgamma(q = 1/maxtausquared, shape = K/2 - 1, rate = 1/2 * crossprod(Theta[g,-c(1, K+2)])) , max = 1)
  # Recovering corresponding inverse tau squared
  invtausquarednew <- qgamma(p = u, shape = K/2 -1, rate = 1/2 * crossprod(Theta[g, -c(1,K+2)]))
  # invert to get new tau squared value
  newtau <- 1/invtausquarednew
  # Update Theta matrix with new tau squared value
  Theta[g,K+2] <- newtau
}

Theta <- Theta[-(1:500),]

plot(Theta[,1], type = "l")
plot(Theta[,2], type = "l")
plot(Theta[,K+2], type = "l")


# Metropolis Hastings Step
# Bring in the Gaussian Approximation with BFGS Optimization
# use log to write the acceptance step
# change index i to g

# Tau Update
# Change update to half-t distribution


