# Create a function of the eta_br model for simulations for single or multiple stationary time series

sim_eta_br = function(n,iter,phi,R=1,B=10){
  
  # Call function to generate data
  data_gen = generate_adapt(phi = phi, n = n, R = R)
  # extract data matrix
  timeseries = data_gen$matrix_timeseries
  
  # call Sampler_eta_r
  return(Sampler_eta_br(timeseries = timeseries, B = 10, iter = 1000, nu = 10, etasq = 1, tausquared = 1, lambda = 1))
}