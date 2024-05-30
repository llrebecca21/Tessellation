#' Function that runs the sampling algorithm for the eta_r model
#'
#' @param timeseries : time series matrix
#' @param B : Number of basis coefficients (excluding the intercept)
#' @param iter : Number of iterations the mcmc algorithm will run
#' @param nu : hyper-parameter for 1/2 t distribution
#' @param etasq : hyper-parameter for 1/2 t distribution
#' @param tausquared : initialization for tau^2
#' @param lambda : initialization for lambda
#'
#' @return
#' @export
#'
#' @examples
Sampler_eta_r = function(timeseries, B = 10, iter = 1000, nu = 3, etasq = 1, tausquared = 1, lambda = 1){
  # Set outer parameters for simulations
  # extract n and R from timeseries
  n = nrow(timeseries)
  R = ncol(timeseries)
  # highest little j index value for the frequencies
  J = floor(n / 2)
  # Frequency (\omega_j): defined on [0, 2\pi)
  omega = (2 * pi * (1:J)) / n
  
  # Define Periodogram
  # Define y_n(\omega_j) for the posterior function below
  perio = (abs(mvfft(timeseries))^2 / n)
  # subset perio for unique values, J = ceil((n-1) / 2) 
  perio = perio[(1:J)+1, , drop=FALSE]
  
  ##########################
  # Set Hyper-Parameters
  ##########################
  # Values that are fixed and set by user
  # Create matrix of the basis functions
  # fix fourier frequencies
  Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
  # redefine the first column to be 1's
  Psi[,1] = 1
  # Define D's main diagonal : 
  # D is a measure of prior variance for \beta_1 through \beta_K
  # Rebecca's D
  D = 1 / (4 * pi * (1:B)^2)
  # Calculate ybar(omega_j)
  y_bar = rowMeans(perio)
  
  ########################
  # Initialize Parameters
  ########################
  # set tau^2 value
  tausquared = 50
  # Define lambda for the half-t prior: \pi(\tau^2 | \lambda) and \pi(\lambda)
  lambda = 1
  # The new D matrix that houses the prior variance of \beta^* 
  Sigma = c(1, D * tausquared)
  # Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
  # ncol: number of parameters (beta^*, tau^2)
  # dim : (iter) x (B + 2)
  Theta = matrix(NA, nrow = iter, ncol = B + 2)
  # Using J amount of data for periodogram, can initialize beta this way:
  betavalues = solve(crossprod(Psi), crossprod(Psi, log(y_bar)))
  # Initialize bb_beta at the mean for the prior of bb_beta
  bb_beta = tcrossprod(betavalues, rep(1,R))
  # Specify Sum of X for the posterior function later
  # Specify Sum of X for the posterior function later
  # 1^T_n X part in the paper: identical to colSums but is a faster calculation
  sumPsi = crossprod(Psi, rep(1,J)) 
  # Initialize first row of Theta
  Theta[1,] = c(betavalues, tausquared)
  # Create Array to store values of bb_beta
  bb_beta_array = array(data = NA, dim = c(iter,nrow(bb_beta),ncol(bb_beta)))
  # initialize first array with bb_beta value
  bb_beta_array[1,,] = bb_beta
  # initialize starting value of eta
  eta = rep(1, length = R)
  # Create Array to store eta values
  eta_array = array(data = NA, dim = c(iter,ncol(bb_beta)))
  # Initialize eta values
  eta_array[1,] = eta
  #####################
  # MCMC Algorithm
  #####################
  # need to update blackboard beta, eta^r, tau^2, lambda 
  for (g in 2:iter) {
    ###############################################
    # Sample \tau^2 and \lambda jointly with Gibbs
    ###############################################
    # Update \lambda and \tau with Gibbs Sampler
    lambda = 1/rgamma(1, (nu+1)/2, nu/tausquared + etasq)
    tausquared = 1/rgamma(1, (nu + B * R)/2, sum(rowSums(bb_beta^2 * tcrossprod(rep(1,B+1),eta))[-1] / D)/2 + nu/lambda)
    # Update Theta matrix with new tau squared value
    Theta[g,B+2] = tausquared
    # Update Sigma with new tau^2 value
    Sigma = c(1, D * tausquared)
    #######################################
    # Sample \eta^r using slicing method
    #######################################
    
    rate = colSums(bb_beta * bb_beta / (2 * Sigma))
    u = runif(R)/(1 + eta)
    q = runif(R)*(pgamma(rate/u^2, shape = (B + 1)/2))
    eta = qgamma(q, shape = (B + 1)/2)/rate
    # put eta_br into array for storage
    eta_array[g,] = eta
    
    ######################
    # bb_beta update :MH
    ######################
    # Update \mathbb{B} with Metropolis Hastings Sampler
    for(r in 1:R){
      #r = 1
      br = bb_beta[,r]
      eta_r = eta[r]
      # Maximum A Posteriori (MAP) estimate : finds the \beta that gives us the mode of the conditional posterior of \beta conditioned on y
      map = optim(par = br, fn = posterior_eta_r, gr = gradient_eta_r, method ="BFGS", control = list(fnscale = -1),
                  Psi = Psi, sumPsi = sumPsi, y = perio[,r], Sigma = Sigma, eta_r = eta_r)$par
      # Call the hessian function
      norm_precision = he_eta_r(br = map, Psi = Psi, y = perio[,r], Sigma = Sigma, eta_r = eta_r) * -2
      # Calculate the \beta^* proposal, using Cholesky Sampling
      betaprop = Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
      # Calculate acceptance ratio
      prop_ratio = min(1, exp(posterior_eta_r(br = betaprop,  Psi = Psi, sumPsi = sumPsi, y = perio[,r], eta_r = eta_r, Sigma = Sigma) -
                                posterior_eta_r(br = br, Psi = Psi, sumPsi = sumPsi, y = perio[,r], eta_r = eta_r, Sigma = Sigma)))
      # Create acceptance decision
      accept <- runif(1)
      if(accept < prop_ratio){
        # Accept betaprop as new beta^(r)
        bb_beta[,r] <- betaprop
      }else{
        # Reject betaprop as new beta^(r)
        bb_beta[,r] <- br
      }
    }
    bb_beta_array[g,,] <- bb_beta
  }
  return(list("bb_beta_array" = bb_beta_array, "eta_array" = eta_array, "Theta" = Theta, "perio" = perio, "av_perio" = y_bar))
  
  
}