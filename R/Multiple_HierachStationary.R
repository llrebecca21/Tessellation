# Multiple Stationary Time Series simulated from perturbations with the same mean from a single time series

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
# create matrix to store the time series: (R x n)
matrix_timeseries = matrix(NA, nrow = n, ncol = R)
# create vector to store "true" theta
theta_vec = rep(NA, R)
for(r in 1:R){
  # set AR parameter
  phi = NULL
  # set MA parameter
  theta = runif(n = 1, min = 0, max = 1)
  theta_vec[r] = theta
  matrix_timeseries[,r] <- arima.sim(model = list(ar = phi, ma = theta), n = n, n.start = burn)
}
dim(matrix_timeseries)
# plot(matrix_timeseries[,1], type = "l")

# Define Periodogram
# Define y_n(\omega_j) for the posterior function below
perio = (abs(mvfft(matrix_timeseries)) ^ 2 / n)
dim(perio)
# subset perio for unique values, J = ceil((n-1) / 2) 
perio = perio[(0:J) + 1, , drop=FALSE]
dim(perio)
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


# \Lambda will be a random piece later
# Define \Lambda for variance of beta^r's
# Lambda = 

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

# Initialize Lambda
Lambda = deg * V

# Specify Sum of X for the posterior function later
# Specify Sum of X for the posterior function later
# 1^T_n X part in the paper: identical to colSums but is a faster calculation
sumPsi = c(crossprod(rep(1, nrow(Psi)), Psi))
# Initialize first row of Theta
Theta[1,] = c(betavalues, tausquared)

# Create Array to store values of bb_beta
bb_beta_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))

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




















