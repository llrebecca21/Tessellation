Sampler_nonstat = function(timeseries, B, iter = 1000){
  # extract n from timeseries
  n = nrow(timeseries)
  # initialization
  # starting break point(s)
  # alpha parameter for Dirichlet distribution
  alpha = rep(1,2)
  # p_vec for Dirichlet prior
  p_vec = rep(.5,2)
  # n_vec for multinomial distribution
  n_vec = c(floor(n/2), ceiling(n/2))
  
  
  for (g in 2:iter) {
    
    
    
  }
  
  
    # call Sampler_eta_br function
    # run the Sampler_eta_br function
    Result_eta_br = Sampler_eta_br(timeseries = timeseries, B = B, tausquared = 1)
    

}