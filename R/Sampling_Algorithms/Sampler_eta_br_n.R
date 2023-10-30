Sampler_eta_br_n = function(ts_list, B = 10, iter = 1000, nu = 3, etasq = 1, tausquared = 1, lambda = 1){
  # extract n and R from timeseries
  # extract the length of each replicate timeseries and store as a vector
  # ts_list will be the input for timeseries
  n_len = sapply(ts_list, nrow)
  R = length(n_len)
  
  # highest little j index value for the frequencies
  J = floor((n_len-1) / 2)
  
  # Define D's main diagonal : 
  # D is a measure of prior variance for \beta_1 through \beta_K
  # Rebecca's D
  D = 1 / (4 * pi * (1:B)^2)
  
  # The new D matrix that houses the prior variance of \beta^* 
  Sigma = c(1, D * tausquared)
  
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
    # Define y_n(\omega_j) for the posterior function below
    perio_list[[r]] = (abs(fft(ts_list[[r]])) ^ 2 / n_len[r])
    
    # subset perio for unique values, J = ceil((n-1) / 2) 
    perio_list[[r]] = perio_list[[r]][(0:J[r]) + 1, , drop = FALSE]
    
    ##########################
    # Set Hyper-Parameter Psi
    ##########################
    # Values that are fixed and set by user
    # Create matrix of the basis functions
    # fix fourier frequencies
    Psi_list[[r]] = outer(X = (2 * pi * (0:J[r])) / n_len[r], Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
    # redefine the first column to be 1's
    Psi_list[[r]][,1] = 1
    
    # Using J amount of data for the periodogram, can initialize beta this way:
    betavalues[,r] = solve(crossprod(Psi_list[[r]]), crossprod(Psi_list[[r]], log(perio_list[[r]])))
    
    # Specify Sum of X for the posterior function later
    # 1^T_n X part in the paper: identical to colSums but is a faster calculation
    #sumPsi[,r] = c(crossprod(rep(1, nrow(Psi_list[[r]])), Psi_list[[r]]))
    sumPsi[,r] = crossprod(Psi_list[[r]], rep(1,J[r]+1)) 
    
  }
  #return(sumPsi)
  betavalues = rowMeans(betavalues)
  
  # Initialize bb_beta at the mean for the prior of bb_beta
  bb_beta = tcrossprod(betavalues, rep(1,R))
  
  # Initialize first row of Theta
  Theta[1,] = c(betavalues, tausquared)
  
  # Create Array to store values of bb_beta
  bb_beta_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
  # initialize first array with bb_beta value
  bb_beta_array[1,,] = bb_beta

  ########################
  # Initialize Parameters
  ########################
  # Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
  # ncol: number of parameters (beta^*, tau^2)
  # dim : (iter) x (B + 2)
  Theta = matrix(NA, nrow = iter, ncol = B + 2)

  # Initialize first row of Theta
  Theta[1,] = c(betavalues, tausquared)
  # Create Array to store values of bb_beta
  bb_beta_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
  # initialize first array with bb_beta value
  bb_beta_array[1,,] = bb_beta
  # initialize starting value of eta_br
  eta_br = matrix(1, nrow = B + 1 , ncol = R)
  # Create Array to store eta_br values
  eta_br_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
  # Initialize eta_br values
  eta_br_array[1,,] = eta_br
  
  #####################
  # MCMC Algorithm
  #####################
  # need to update blackboard beta, eta_b^r, tau^2, lambda 
  for (g in 2:iter) {
    ###############################################
    # Sample \tau^2 and \lambda jointly with Gibbs
    ###############################################
    # Update \lambda and \tau with Gibbs Sampler
    lambda = 1/rgamma(1, (nu+1)/2, nu/tausquared + etasq)
    tausquared = 1/rgamma(1, (nu + B * R)/2, sum(rowSums(bb_beta^2 * eta_br)[-1] / D)/2 + nu/lambda)
    # Update Theta matrix with new tau squared value
    Theta[g,B+2] = tausquared
    # Update Sigma with new tau^2 value
    Sigma = c(1, D * tausquared)
    #######################################
    # Sample \eta_b^r using slicing method
    #######################################
    rate = c(bb_beta * bb_beta / (2 * Sigma))
    p = R * (B + 1)
    u = runif(p)/(1 + c(eta_br))
    q = runif(p)*(pgamma(rate/u^2, shape = 1))
    eta_br = qgamma(q, shape = 1)/rate
    # reconstruct shape
    eta_br = matrix(eta_br, nrow = B + 1 , ncol = R)
    # put eta_br into array for storage
    eta_br_array[g,,] = eta_br
    
    ######################
    # bb_beta update :MH
    ######################
    # Update \mathbb{B} with Metropolis Hastings Sampler
    for(r in 1:R){
      #r = 1
      br = bb_beta[,r]
      eta_r = eta_br[,r]
      # Maximum A Posteriori (MAP) estimate : finds the \beta that gives us the mode of the conditional posterior of \beta conditioned on y
      map = optim(par = br, fn = posterior_eta_br, gr = gradient_eta_br, method ="BFGS", control = list(fnscale = -1),
                  Psi = Psi_list[[r]], sumPsi = sumPsi[ , r, drop = FALSE], y = perio_list[[r]], Sigma = Sigma, eta_r = eta_r)$par
      # Call the hessian function
      norm_precision = he_eta_br(br = map, Psi = Psi_list[[r]], y = perio_list[[r]], Sigma = Sigma, eta_r = eta_r) * -1
      # Calculate the \beta^* proposal, using Cholesky Sampling
      betaprop = Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
      # Calculate acceptance ratio
      prop_ratio = min(1, exp(posterior_eta_br(br = betaprop,  Psi = Psi_list[[r]], sumPsi = sumPsi[,r,drop = FALSE], y = perio_list[[r]], eta_r = eta_r, Sigma = Sigma) -
                                posterior_eta_br(br = br, Psi = Psi_list[[r]], sumPsi = sumPsi[,r,drop = FALSE], y = perio_list[[r]], eta_r = eta_r, Sigma = Sigma)))
      # Create acceptance decision
      accept <- runif(1)
      if(accept < prop_ratio){
        # Accept betaprop as new beta^(r)
        bb_beta_array[g, ,r] <- betaprop
      }else{
        # Reject betaprop as new beta^(r)
        bb_beta_array[g, ,r] <- br
      }
    }
    
  }
  return(list("bb_beta_array" = bb_beta_array, "eta_br_array" = eta_br_array, "Theta" = Theta, "perio_list" = perio_list))
}