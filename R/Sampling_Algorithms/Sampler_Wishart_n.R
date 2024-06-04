Sampler_Wishart_n = function(ts_list,  B = 10, iter = 1000, nu = 3, etasq = 1, sigmasquare = 100, tausquared = 1, lambda = 1, V = NULL){
  
  # extract n and R from timeseries
  # extract the length of each replicate timeseries and store as a vector
  # ts_list will be the input for timeseries
  n_len = sapply(ts_list, nrow)
  R = length(n_len)
  
  # highest little j index value for the frequencies
  J = floor(n_len / 2)
  # Frequency (\omega_j): defined on [0, 2\pi)
  #omega = (2 * pi * (0:J)) / n_len
  
  #################
  # MCMC parameters
  #################
  # Define D's main diagonal : 
  # D is a measure of prior variance for \beta_1 through \beta_K
  # Rebecca's D
  D = 1 / (4 * pi * (1:B)^2)
  
  # prior for Lambda : Wishart
  # degrees of freedom
  deg = B+1
  if(is.null(V)){
    # initial V
    V = diag(B+1)
  }
  
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
    perio_list[[r]] = perio_list[[r]][(1:J[r]) + 1, , drop = FALSE]
    
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
    #sumPsi[,r] = c(crossprod(rep(1, nrow(Psi_list[[r]])), Psi_list[[r]]))
    sumPsi[,r] = crossprod(Psi_list[[r]], rep(1,J[r])) 
  }
  #return(sumPsi)
  betavalues = rowMeans(betavalues)
  
  # Initialize bb_beta at the mean for the prior of bb_beta
  bb_beta = tcrossprod(betavalues, rep(1,R))
  
  # Initialize first row of Theta
  Theta[1,] = c(betavalues, tausquared)
  
  # Initialize Lambda^{-1}
  Lambda_inv = deg * V
  
  # Create Array to store values of bb_beta
  bb_beta_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
  # initialize first array with bb_beta value
  bb_beta_array[1,,] = bb_beta
  
  # Create array to hold Lambda values
  Lambda_array = array(data = NA, dim = c(iter, nrow(Lambda_inv), ncol(Lambda_inv)))
  # Initalize first array with Lambda value
  Lambda_array[1,,] = Lambda_inv
  
  #####################
  # MCMC Algorithm
  #####################
  
  for (g in 2:iter) {
    #g = 2
    #########################
    # tau^2 and lambda update
    #########################
    # Update \lambda and \tau with Gibbs Sampler
    lambda = 1/rgamma(1, (nu+1)/2, nu/tausquared + etasq)
    tausquared = 1/rgamma(1, (B + 1 + nu)/2, sum(betavalues[-1]^2 / D) / 2 + nu/lambda)
    # Update Theta matrix with new tau squared value
    Theta[g,B+2] = tausquared
    # Update Sigma with new tau^2 value
    Sigma = c(sigmasquare, D * tausquared)
    
    #####################
    # beta update
    #####################
    # Update \beta with Gibbs Sampler
    Lambda_solve = solve(diag(Sigma))
    betavalues = c(mvtnorm::rmvnorm(n = 1, mean = solve(Lambda_solve + R * Lambda_inv) %*% Lambda_inv %*% bb_beta %*% rep(1,R),
                           sigma = solve(Lambda_solve + R * Lambda_inv)))
    #save new betavalue
    Theta[g, -(B+2)] = betavalues
    
    #####################
    # Lambda^{-1} update
    #####################
    # Update \Lambda^{-1} with Gibbs Sampler
    # obtain new Lambda_inv with Wishart distribution
    Lambda_inv = rWishart(n = 1, df = deg + R, Sigma = solve(tcrossprod(bb_beta - tcrossprod(betavalues, rep(1,R))) + diag(B+1)))[,,1]
    # save Lambda_inv 
    Lambda_array[g,,] = Lambda_inv
    
    ######################
    # bb_beta update :MH
    ######################
    # Update \mathbb{B} with Metropolis Hastings Sampler
    for(r in 1:R){
      #r = 1
      br = bb_beta[,r]
      # Maximum A Posteriori (MAP) estimate : finds the \beta that gives us the mode of the conditional posterior of \beta conditioned on y
      map <- optim(par = br, fn = posterior_hierarch_Lambda, gr = gradient_hierarch_Lambda, method ="BFGS", control = list(fnscale = -1),
                   Psi = Psi_list[[r]], sumPsi = sumPsi[,r, drop = FALSE], y = perio_list[[r]], b = betavalues, lambda = Lambda_inv)$par
      # Call the hessian function
      norm_precision <- he_hierarch_Lambda(br = map, Psi = Psi_list[[r]], y = perio_list[[r]], lambda = Lambda_inv) * -1
      # Calculate the \beta^* proposal, using Cholesky Sampling
      betaprop <- Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
      # Calculate acceptance ratio
      prop_ratio <- min(1, exp(posterior_hierarch_Lambda(br = betaprop, b = betavalues, Psi = Psi_list[[r]], sumPsi = sumPsi[,r, drop = FALSE], y = perio_list[[r]], lambda = Lambda_inv) -
                                 posterior_hierarch_Lambda(br = br, b = betavalues, Psi = Psi_list[[r]], sumPsi = sumPsi[,r, drop = FALSE], y = perio_list[[r]], lambda = Lambda_inv)))
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
  
  return(list("bb_beta_array" = bb_beta_array, "Lambda_array" = Lambda_array, "Theta" = Theta, "perio_list" = perio_list))
  
}