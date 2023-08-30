# Multiple Stationary Time Series simulated from perturbations with the same mean from a single time series

# Create time series 
# set hyper-parameters
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





















