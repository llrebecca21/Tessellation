# Create a function of the eta_r model for simulations for single or multiple stationary time series

sim_eta_r = function(n,iter,phi,R,B,nu,etasq,tausquared,lambda,theta_z,nu_r,burnin, scaling){
  # check the correct values are fed into the Sampler

  # Call function to generate data
  data_gen = generate_adapt(phi = phi, n = n, R = R)
  # extract data matrix
  timeseries = data_gen$matrix_timeseries
  
  # call Sampler_eta_r
  return(Sampler_eta_r(timeseries, B, iter, nu = 10, etasq = 1, tausquared = 1, lambda = 1, theta_z = 1, nu_r = 1, burnin = burnin, scaling = 1))
}

















