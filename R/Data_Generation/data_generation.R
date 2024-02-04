generate_Krafty = function(n = 1000, R = 1,burn = 50){
  # create matrix to store the time series: (R x n)
  matrix_timeseries = matrix(NA, nrow = n, ncol = R)
  # create vector to store "true" theta
  theta_true = rep(NA, R)
  for(r in 1:R){
    # set AR parameter
    phi = NULL
    # set MA parameter
    theta = runif(n = 1, min = 0, max = 1)
    theta_true[r] = theta
    matrix_timeseries[,r] <- arima.sim(model = list(ar = phi, ma = theta), n = n, n.start = burn)
  }
  return(list("matrix_timeseries" = matrix_timeseries, "theta_true" = theta_true))
}

# Generate a general AR(p) process
generate_adapt = function(phi, n = 1000, R = 1, burn = 50){
  # Need to Create ~ R copies of the time series and store it in a matrix
  # Each column of the matrix contains a time series
  # create matrix to store the time series: (R x n)
  matrix_timeseries = matrix(NA, nrow = n, ncol = R)
  for(r in 1:R){
    matrix_timeseries[,r] <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)
  }
  return(list("matrix_timeseries" = matrix_timeseries))
}


# Generate an AR(2) with angle representation
generate_ar2_peak = function(peaks, bandwidths, variances = NULL, n){
  # peaks - vector of length R with peak locations of each spectrum
  # bandwidths - vector of length n_ts with the bandwidths of each spectrum
  # variances - vector of length n_ts with the variances of the innovations 
  #   of each process
  # n_pts - length of each realization to be drawn of the process
  R <- length(peaks)
  if(is.null(variances)){
    variances = rep(1, R)
  }
  # generate matrix of realizations of each process, the multitaper estimates,
  # and the "true" spectra
  matrix_timeseries <- matrix(nrow = n, ncol = R)
  phi <- matrix(NA, nrow = 2, ncol = R)
  for(i in 1:R) {
    # autoregressive parameters
    phi[1,i] <- 2 * cos(peaks[i]) * exp(-bandwidths[i])
    phi[2,i] <- -exp(-2 * bandwidths[i])
    
    # data
    matrix_timeseries[, i] <- arima.sim(list(ar = phi[,i]), 
                        n = n, 
                        sd = sqrt(variances[i]))
    
  }
  return(list(
    "matrix_timeseries" = matrix_timeseries,
    "phi" = phi
  ))
}


# Generate a non-stationary time series
generate_nonstat_abrupt = function(phi1, phi2, n = 1000, R = 1, burn = 50){
  # phi1 controls the input for the first AR(p) process
  # phi2 controls the input for the second AR(p) process
  # bp : controls how many breakpoints we want
  
  # Need to Create ~ R copies of the time series and store it in a matrix
  # Each column of the matrix contains a time series
  # create matrix to store the time series: (R x (nx2))
  matrix_timeseries = matrix(NA, nrow = n * bp, ncol = R)
  for(r in 1:R){
    matrix_timeseries[,r] <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)
  }
  return(list("matrix_timeseries" = matrix_timeseries))
}








