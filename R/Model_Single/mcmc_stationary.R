#' Function that rungs the mcmc sampling algorithm for a single or multiple stationary time series
#'
#' @param n : (integer) length of an individual time series
#' @param iter : (integer) number of iterations the mcmc algorithm completes
#' @param R : (integer) number of time series from same underlying spectrum
#' @param B : (integer) number of basis functions
#' @param phi : (vector or value between 0 and 1) coefficient value(s) of the AR process used to create data
#'
#' @return
#' @export
#'
#' @examples
mcmc_stationary <- function(n, iter, phi ,R = 1, B = 10 ){
  
  # Create a single time series
  # set hyper-parameters
  # highest little j index value for the frequencies
  J = floor((n-1) / 2)
  # Frequency (\omega_j): defined on [0, 2\pi)
  # for j = 0,...,n-1
  # omega = (2 * pi * (0:(n-1)))/n
  omega = (2 * pi * (0:J)) / n
  
  # Call function to generate data
  data_gen = generate_adapt(phi = phi, n = n, R = R)
  # extract data matrix
  matrix_timeseries = data_gen$matrix_timeseries

  # Define Periodogram
  # Define y_n(\omega_j) for the posterior function below
  perio = (abs(mvfft(matrix_timeseries)) ^ 2 / n)

  # subset perio for unique values, J = ceil((n-1) / 2) 
  perio = perio[(0:J) + 1, , drop=FALSE]

  # Calculate ybar(omega_j)
  y_bar = rowMeans(perio)
  
  #################
  # MCMC parameters
  #################

  # Define lambda for the half-t prior: \pi(\tau^2 | \lambda) and \pi(\lambda)
  lambda = 1
  # Define degrees of freedom for half-t prior: \pi(\tau^2 | \lambda)
  # Cauchy nu = 1
  nu = 3
  # Define eta as the other scale parameter for half-t prior: \pi(\lambda)
  # etasq = 1 gives standard Cauchy; higher eta gives wider Cauchy
  etasq = 1
  
  # Define D's main diagonal : Choose either identity or Rebecca's D (variation of Yakun's D)
  # D is a measure of prior variance for \beta_1 through \beta_K
  
  # Identity D:
  # D = rep(1, B)
  
  # Rebecca's D
  D = 1 / (4 * pi * (1:B)^2)
  
  # exponential decay D
  # D = exp(0.12 * -(0:(B-1)))
  
  # prior variance for beta_0
  sigmasquare = 100
  
  #######################
  # Initialize parameters
  #######################
  # set tau^2 value
  tausquared = 1
  # The new D matrix that houses the prior variance of \beta^* 
  Sigma = c(sigmasquare, D * tausquared)
  
  # Create matrix to store estimated samples row-wise for (\beta^*, \tau^2)
  # ncol: number of parameters (beta^*, tau^2)
  # dim : (iter) x (B + 2)
  Theta = matrix(NA, nrow = iter, ncol = B + 2)
  
  # Create matrix of the basis functions
  # fix fourier frequencies
  Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
  # redefine the first column to be 1's
  Psi[,1] = 1
  # not orthogonal because we are not evaluating the periodogram at the full n-1 values.
  # Initialize beta using least squares solution
  # Using the full data n-1 for periodogram can use this to initialize beta:
  #betavalues = c(crossprod(Psi,log(y_bar))) / n
  
  # Using J amount of data for periodogram, can initialize beta this way:
  betavalues = solve(crossprod(Psi), crossprod(Psi, log(y_bar)))
  
  
  # Specify Sum of X for the posterior function later
  # Specify Sum of X for the posterior function later
  # 1^T_n X part in the paper: identical to colSums but is a faster calculation
  sumPsi = c(crossprod(rep(1, nrow(Psi)), Psi))
  # Initialize first row of Theta
  Theta[1,] = c(betavalues, tausquared)
  
  #####################
  # MCMC Algorithm
  #####################
  
  #Rprof()
  # pb = progress_bar$new(total = iter - 1)
  for (g in 2:iter) {
    # pb$tick()
    # g = 2
    # Extract \beta^* and tau^2 from theta
    # beta^* of most recent iteration:
    b = Theta[g - 1, 1:(B+1)]
    # tau^squared of most recent iteration:
    tsq = Theta[g - 1, B + 2]
    ##########################
    # Metropolis Hastings Step
    ##########################
    # Maximum A Posteriori (MAP) estimate : finds the \beta^* that gives us the mode of the conditional posterior of \beta^* conditioned on y
    map <- optim(par = b, fn = posterior_multiple, gr = gr_multiple, method ="BFGS", control = list(fnscale = -1),
                 Psi = Psi, sumPsi = sumPsi, y_bar = y_bar, D = Sigma, R = R)$par
    # Call the hessian function
    norm_precision <- he_multiple(b = map, Psi = Psi, y_bar = y_bar, D = Sigma , R = R)
    # Calculate the \beta^* proposal, using Cholesky Sampling
    betaprop <- Chol_sampling(Lt = chol(norm_precision), d = B + 1, beta_c = map)
    # Calculate acceptance ratio
    prop_ratio <- min(1, exp(posterior_multiple(b = betaprop, Psi = Psi, sumPsi = sumPsi, y_bar = y_bar,  D = Sigma, R = R) -
                               posterior_multiple(b = b, Psi = Psi, sumPsi = sumPsi, y_bar = y_bar,  D = Sigma, R = R)))
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
  return(list("Theta" = Theta, "av_perio" = y_bar))
  
  
  

}






