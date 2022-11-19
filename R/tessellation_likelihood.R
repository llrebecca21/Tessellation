### Whittle likelihood for each tessellation
# y is a list of log periodogram for ith tessellation and lth time series with T[i,l] length, T[i,l] is the ith tessellation for lth time series
# y is y[i][[l]]
# fhat is a list of log power spectrum for ith tessellation and lth time series with T[i,l] length
# beta is the matrix of collection of basis parameters for each tessellation beta[i,j]: ith tessellation, jth basis 

#### whittle_like is for ith tessellation: T is T[i,] and beta[i,]

whittle_like <-
  function(y,index,nseg_time_temp,beta,nbasis) {
    
    log_prop_spec_dens=0
    nbeta=nbasis+1
    
    uu=0
    for(l in index)
    {
      
      uu=uu+1
      nfreq=floor(nseg_time_temp[l]/2)
      freq <- (0:nfreq)/(2 * nfreq)
      nu_mat <- lin_basis_func(freq, nbeta)
      fhat <- nu_mat%*%matrix(beta,nbeta,1)
      
      
      if(nseg_time_temp[l]%%2==1)
      {
        
        log_prop_spec_dens=log_prop_spec_dens-sum(fhat[2:(nfreq+1)]+exp(y[[uu]][2:(nfreq+1)]-fhat[2:(nfreq+1)]))-
          0.5*(fhat[1]+exp(y[[uu]][1]-fhat[1]))-0.5*nfreq*log(2*pi)
        
      }else
      {
        log_prop_spec_dens=log_prop_spec_dens-sum(fhat[2:nfreq]+exp(y[[uu]][2:nfreq]-fhat[2:nfreq]))-
          0.5*(fhat[1]+exp(y[[uu]][1]-fhat[1]))-0.5*(fhat[nfreq+1]+exp(y[[uu]][nfreq+1]-fhat[nfreq+1]))-
          0.5*nfreq*log(2*pi)
      }
    }
    
    return(log_prop_spec_dens)
    
    
  }



# whittle_derivs2 is for the ith tessellation
# param is the initialize of beta=(0,0,0,0)
# y[[uu]] is log periodogram for the lth time series
# 
whittle_derivs2 <-
  function(param,
           index,
           nseg_time_temp,
           y,
           tau_temp,
           nbeta,
           nbasis,
           sigmasqalpha)
  {
    uu=0
    for(l in index)
    {
      uu=uu+1
      nfreq <- floor(nseg_time_temp[l]/2)
      freq <- (0:nfreq)/(2 * nfreq)
      nu_mat <- lin_basis_func(freq, nbeta)
      ydev <- y[[uu]]-nu_mat %*% param     ### ydev is a vector
      
      n1 <- nfreq
      
      f=0
      
      if (nseg_time_temp[l] %% 2 == 1)
        # odd n
      {
        
        f <- f + sum(nu_mat[2:(n1 + 1), ] %*% param + exp(ydev[2:(n1 + 1)])) +
          0.5 * (sum(as.vector(nu_mat[1, ] %*% param) + exp(ydev[1]))) +
          0.5 * (t(param[2:nbeta]) %*% param[2:nbeta] / tau_temp + param[1] ^
                   2 / sigmasqalpha)
        
      } else
      {
        f <- f + sum(nu_mat[2:n1, ] %*% param + exp(ydev[2:n1])) +
          0.5 * (sum(as.vector(nu_mat[1, ] %*% param) + exp(ydev[1]))) +
          0.5 * (sum(as.vector(nu_mat[n1 + 1, ] %*% param) + exp(ydev[n1 + 1]))) +
          0.5 * (t(param[2:nbeta]) %*% param[2:nbeta] / tau_temp + param[1] ^
                   2 / sigmasqalpha)
      }
      
    }
    
    g <- matrix(0, nbeta, 1)
    g[1, ] <- param[1] / sigmasqalpha
    g[2:nbeta, ] <- param[2:nbeta] / tau_temp
    
    h <- matrix(0, nbeta, nbeta)
    h[1, 1] <- 1 / sigmasqalpha
    h[2:nbeta, 2:nbeta] <- 1 / tau_temp * diag(nbasis)
    
    
    uu=0
    for (l in index)
    {
      uu=uu+1 
      if (nseg_time_temp[l] %% 2 == 1)
      {
        
        
        nfreq <- floor(nseg_time_temp[l]/2)
        freq <- (0:nfreq)/(2 * nfreq)
        nu_mat <- lin_basis_func(freq, nbeta)
        ydev <- y[[uu]]-nu_mat %*% param     ### ydev is a vector
        
        n1 <- nfreq
        
        temp_mat <- nu_mat[2:(n1 + 1), ] * repmat((1 - exp(as.matrix(ydev[2:(n1 + 1)]))), 1, dim(nu_mat)[2])
        g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1])) %*%
                                                          nu_mat[1, ])
        
        
        jj <- seq(1:nbeta)
        big_mat <- repmat(t(nu_mat[2:(n1 + 1), ]), nbeta, 1) * t(nu_mat[2:(n1 + 1), repmat(jj, nbeta, 1)])
        
        coef_mat <- repmat(exp(t(ydev[2:(n1 + 1)])), nbeta ^ 2, 1)
        h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1)), c(1, 2), sum) +
          t(0.5 * exp(ydev[1]) %*% nu_mat[1, ]) %*% nu_mat[1, ]
        
        
        
      }else{
        
        nfreq <- floor(nseg_time_temp[l]/2)
        freq <- (0:nfreq)/(2 * nfreq)
        nu_mat <- lin_basis_func(freq, nbeta)
        ydev <- y[[uu]]-nu_mat %*% param     ### ydev is a vector
        
        n1 <- nfreq
        
        temp_mat <- nu_mat[2:n1, ] * repmat((1 - exp(as.matrix(ydev[2:n1]))), 1, dim(nu_mat)[2])
        g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1])) %*% nu_mat[1, ]) +
          t(0.5 * (1 - exp(ydev[n1 + 1])) %*% nu_mat[n1 + 1, ])
        
        
        jj <- seq(1:nbeta)
        big_mat <- repmat(t(nu_mat[2:n1, ]), nbeta, 1) * t(nu_mat[2:n1, repmat(jj, nbeta, 1)])
        
        
        coef_mat <- repmat(exp(t(ydev[2:n1])), nbeta ^ 2, 1)
        h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1 - 1)), c(1, 2), sum) +
          t(0.5 * exp(ydev[1]) %*% nu_mat[1, ]) %*% nu_mat[1, ] +
          t(0.5 * exp(ydev[n1 + 1]) %*% nu_mat[n1 + 1, ]) %*% nu_mat[n1 + 1, ]
        
      }
      
    }
    
    
    list(value = f,
         gradient = g,
         hessian = h)
    
  }



lin_basis_func <-
  function(freq, nbeta)
  {
    n <- length(freq)
    omega <- matrix(1, n, nbeta)
    for (j in 2:nbeta)
    {
      omega[, j] <- sqrt(2) * cos((j-1) * pi * freq)/(pi*(j-1))
    }
    
    return(omega)
  }


# i is for the ith tessellation
# index is for the index of time series in the ith tessellation

postbeta <- function(i,index,tmax,tmin,nseg_time_temp,x_t,tau_temp,nbasis,sigmasqalpha)
{
  
  
  # create log periodograms for replications in ith time seg and jth cov seg
  y <- list()
  yy <- list()
  
  uu=0
  for(l in index)
  {
    uu=uu+1
    
    nfreq <- floor(nseg_time_temp[l] / 2)
    
    y[[uu]] <- log(abs(fft(x_t[tmin[l]:tmax[l], l])) ^ 2 / nseg_time_temp[l])
    yy[[uu]] <- y[[uu]][1:(nfreq + 1)]
    
  }
  
  # pass log periodograms into optimizer to obtain beta_mean and beta_var for normal approximation to beta posterior
  nbeta <- nbasis + 1
  # nu_mat <- lin_basis_func(freq, nbeta) # basis function
  param <- rep(0, nbeta)
  
  post <- trust(whittle_derivs2, param, rinit = 1, rmax = 100,
                parscale = rep(1, nbeta), iterlim = 100, fterm = sqrt(.Machine$double.eps),
                mterm = sqrt(.Machine$double.eps), minimize = TRUE,blather = FALSE,
                index, nseg_time_temp, yy, tau_temp, nbeta, nbasis,sigmasqalpha)
  
  beta_mean <- post$argument
  beta_var <-  solve(post$hessian)
  list(beta_mean = beta_mean, beta_var = beta_var, nu_mat = nu_mat,yy=yy,
       y = y)
  
}