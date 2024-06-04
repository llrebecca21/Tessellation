Sampler_Single_n = function(ts_list, B = 10, iter = 1000, nu = 3, etasq = 1, sigmasquare = 100, tausquared = 1, lambda = 1){

  
  # extract n and R from timeseries
  # extract the length of each replicate timeseries and store as a vector
  # ts_list will be the input for timeseries
  n_len = sapply(ts_list, nrow)
  R = length(n_len)
  
  # Create a single time series
  # set hyper-parameters
  # highest little j index value for the frequencies
  J = floor(n_len / 2)
  # Frequency (\omega_j): defined on [0, 2\pi)
  # for j = 0,...,n-1
  # omega = (2 * pi * (0:J)) / n
  
  #################
  # MCMC parameters
  #################
  
  # Rebecca's D
  D = 1 / (4 * pi * (1:B)^2)
  
  #######################
  # Initialize parameters
  #######################
  # The new D matrix that houses the prior variance of \beta^* 
  Sigma = c(sigmasquare, D * tausquared)
  
  # Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
  # ncol: number of parameters (beta^*, tau^2)
  # dim : (iter) x (B + 2)
  Theta = matrix(NA, nrow = iter, ncol = B + 2)
  # initialize perio as a list
  perio_list = vector(mode = "list", length = R)
  # initialize Psi as a list
  Psi_list = vector(mode = "list", length = R)
  # initialize beta_values
  betavalues = matrix(data = NA, nrow = B + 1, ncol = R)
  # initialize sumPsi
  sumPsi = matrix(data = NA, nrow = B + 1, ncol = R)
  
  # Define Periodogram and Psi
  for(r in 1:R){
    #r = 1
    # Define y_n(\omega_j) for the posterior function below
    perio_list[[r]] = (abs(fft(ts_list[[r]]))^ 2 / n_len[r])
    
    # subset perio for unique values, J = ceil((n-1) / 2) 
    perio_list[[r]] = perio_list[[r]][(1:J[r]+1), , drop = FALSE]
    
    # Create matrix of the basis functions
    # fix fourier frequencies
    Psi_list[[r]] = outer(X = (2 * pi * (1:J[r])) / n_len[r], Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
    # redefine the first column to be 1's
    Psi_list[[r]][,1] = 1
    
    # Using J amount of data for periodogram, can initialize beta this way:
    betavalues[,r] = solve(crossprod(Psi_list[[r]]), crossprod(Psi_list[[r]], log(perio_list[[r]])))
    
    # Specify Sum of X for the posterior function later
    # Specify Sum of X for the posterior function later
    # 1^T_n X part in the paper: identical to colSums but is a faster calculation
    sumPsi[,r] = c(crossprod(rep(1, nrow(Psi_list[[r]])), Psi_list[[r]]))
  }
  betavalues = rowMeans(betavalues)
  # calculate sumsumPsi which is the row sums of sumPsi
  sumsumPsi = rowSums(sumPsi)

  # Initialize first row of Theta
  Theta[1,] = c(betavalues, tausquared)
  
  #####################
  # MCMC Algorithm
  #####################
  
  #Rprof()
  # pb = progress_bar$new(total = iter - 1)
  for (g in 2:iter) {
      # pb$tick()
      #g = 2
      # Extract \beta^* and tau^2 from theta
      # beta^* of most recent iteration:
      b = Theta[g - 1, 1:(B+1)]
      # tau^squared of most recent iteration:
      tsq = Theta[g - 1, B + 2]
      ##########################
      # Metropolis Hastings Step
      ##########################
      # Maximum A Posteriori (MAP) estimate : finds the \beta^* that gives us the mode of the conditional posterior of \beta^* conditioned on y
      map = optim(par = b, fn = posterior_multiple_n, gr = gr_multiple_n, method ="BFGS", control = list(fnscale = -1),
                   Psi_list = Psi_list, sumsumPsi = sumsumPsi, perio_list = perio_list, Sigma = Sigma, R = R)$par
      # Call the hessian function, and multiple by -1 to make it a positive definite matrix
      norm_precision <- -1 * he_multiple_n(b = map, Psi_list = Psi_list, perio_list = perio_list, Sigma = Sigma , R = R)
      # Calculate the \beta^* proposal, using Cholesky Sampling
      betaprop <- Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
      # Calculate acceptance ratio
      prop_ratio <- min(1, exp(posterior_multiple_n(b = betaprop, Psi_list = Psi_list, sumsumPsi = sumsumPsi, perio_list = perio_list,  Sigma = Sigma, R = R) -
                                 posterior_multiple_n(b = b, Psi_list = Psi_list, sumsumPsi = sumsumPsi, perio_list = perio_list,  Sigma = Sigma, R = R)))
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
    
    
    #######################
    # Plots and Diagnostics
    #######################
    # Remove burn-in
    burnin <- 10
    Theta <- Theta[-(1:burnin),]
    return(list("Theta" = Theta, "perio" = perio_list))
    
  
}