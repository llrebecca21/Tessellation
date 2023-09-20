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