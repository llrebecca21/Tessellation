Sampler_nonstat = function(timeseries, B, iter = 1000, nmin = 30){
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
  iter_beta = 5
  
  
  for (g in 2:iter) {
    
    # call Whittle post
    # store values for p in a vector (transition probabilities)
    P =  rep(0, n+1)
    
    for(n1 in 0:n){
      n2 = n-n1
      # Get updated betas from Sampler_eta_br
      # Result_eta_br here returns a vector
      out_1 = Sampler_eta_br(timeseries = timeseries[1:n1, ,drop=FALSE], B = B, tausquared = 1, iter = iter_beta)
      out_2 = Sampler_eta_br(timeseries = timeseries[(n1+1):n, ,drop=FALSE], B = B, tausquared = 1, iter = iter_beta)
      #extract what we need from the outputs
      beta_1 = out_1$bb_beta_array[iter_beta,,1]
      beta_2 = out_2$bb_beta_array[iter_beta,,1]
      Sigma_1 = out_1$Sigma
      Sigma_2 = out_2$Sigma
      Psi_1 = out_1$Psi
      Psi_2 = out_2$Psi
      sumPsi_1 = out_1$sumPsi
      sumPsi_2 = out_2$sumPsi
      eta_r_1 = out_1$eta_br_array[iter_beta, ,1]
      eta_r_2 = out_2$eta_br_array[iter_beta, ,1]
      perio_1 = out_1$perio[,1]
      perio_2 = out_2$perio[,1]
      
      
      # Whittle Posterior
      # highest little j index value for the frequencies
      J = floor((n1-1) / 2)
      # Frequency (\omega_j): defined on [0, 2\pi)
      omega = (2 * pi * (0:J)) / n1
      # Create Psi_1
      # fix fourier frequencies
      Psi_1 = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
      # redefine the first column to be 1's
      Psi_1[,1] = 1
      
      # Create Psi_2
      # highest little j index value for the frequencies
      J = floor((n2-1) / 2)
      # Frequency (\omega_j): defined on [0, 2\pi)
      omega = (2 * pi * (0:J)) / n2
      # fix fourier frequencies
      Psi_2 = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
      # redefine the first column to be 1's
      Psi_2[,1] = 1
      
      
      
    }
    
  }
  
  

    

}