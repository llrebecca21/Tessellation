# Create a function of the independent model for simulations for single or multiple stationary time series

sim_independent = function(n,iter,phi,R=1,B=10){
  
  # Call function to generate data
  data_gen = generate_adapt(phi = phi, n = n, R = R)
  # extract data matrix
  timeseries = data_gen$matrix_timeseries
  
  # call Sampler_eta_r
  return(Sampler_Single(timeseries = timeseries, B = 10, iter = 1000, nu = 10, etasq = 1, sigmasquare = 100, tausquared = 1, lambda = 1, burnin = 50))
}