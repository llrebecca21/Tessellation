fn=function(param,
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
  
  return(f)
  
}



gr=function(param,
            index,
            nseg_time_temp,
            y,
            tau_temp,
            nbeta,
            nbasis,
            sigmasqalpha)
{
  
  
  g <- matrix(0, nbeta, 1)
  g[1, ] <- param[1] / sigmasqalpha
  g[2:nbeta, ] <- param[2:nbeta] / tau_temp
  
  
  
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
      
      
    }else{
      
      nfreq <- floor(nseg_time_temp[l]/2)
      freq <- (0:nfreq)/(2 * nfreq)
      nu_mat <- lin_basis_func(freq, nbeta)
      ydev <- y[[uu]]-nu_mat %*% param     ### ydev is a vector
      
      n1 <- nfreq
      
      temp_mat <- nu_mat[2:n1, ] * repmat((1 - exp(as.matrix(ydev[2:n1]))), 1, dim(nu_mat)[2])
      g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1])) %*% nu_mat[1, ]) +
        t(0.5 * (1 - exp(ydev[n1 + 1])) %*% nu_mat[n1 + 1, ])
      
      
    }
    
  }
  
  return(g)
  
  
}


he=function(param,
            index,
            nseg_time_temp,
            y,
            tau_temp,
            nbeta,
            nbasis,
            sigmasqalpha)
{
  
  
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
      
      
      jj <- seq(1:nbeta)
      big_mat <- repmat(t(nu_mat[2:n1, ]), nbeta, 1) * t(nu_mat[2:n1, repmat(jj, nbeta, 1)])
      
      
      coef_mat <- repmat(exp(t(ydev[2:n1])), nbeta ^ 2, 1)
      h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1 - 1)), c(1, 2), sum) +
        t(0.5 * exp(ydev[1]) %*% nu_mat[1, ]) %*% nu_mat[1, ] +
        t(0.5 * exp(ydev[n1 + 1]) %*% nu_mat[n1 + 1, ]) %*% nu_mat[n1 + 1, ]
      
    }
    
  }
  
  return(h)
  
}


###### Cholesky sampling
Chol_sampling=function(Q,d,mean)
{
  L=t(chol(Q))   ## L : left Cholesky factor
  Lt=chol(Q)    ## right Cholesky factor
  z=t(rmvnorm(1,rep(0,d),diag(d)))
  u=solve(Lt,z)
  beta_prop=mean+u
  
  return(beta_prop)
  
}


###### Cholesky density
Chol_density=function(x,Q,d,mean)
{
  
  L=t(chol(Q))
  Lt=chol(Q)
  x_tilde=Lt%*%(x-mean)
  term1=(-d/2)*log(2*pi)
  term2=sum(log(diag(L)))
  term3=(-1/2)*sum(x_tilde^2)
  log_beta_prop=term1+term2+term3
  
  return(log_beta_prop)
  
}


## Whithin step
Whithin=function(x,x_t,TT,NumObs,nbasis,sigmasqalpha,pre_beta)
{
  
  p=dim(x)[2]-2  # remove index and time
  n=dim(x)[1]    # number of observations
  nbeta=nbasis+1
  
  
  # proposed partition
  prt=distance_partitionC(as.matrix(x[,-c(1,2)]),TT$S,TT$w)
  
  # time interval
  interval_curr=time_interval(x,TT$S,prt,NumObs)
  
  
  
  #### proposed model
  
  # Drawing new tausq and g
  tau_prop=TT$tau
  g_prop=TT$g
  
  
  ## calculate log power spectrum
  nfreq=50  # fixed the number of Fourier frequencies
  freq <- (0:nfreq)/(2 * nfreq)
  nu_mat <- lin_basis_func(freq, nbeta)
  
  
  
  alpha=rep(0,TT$M)
  accepted=rep(0,TT$M)
  
  # update beta and loglike
  for(i in 1:TT$M)
  {
    
    
    prop_obs=x[which(prt==i),]
    prop_obs_index=unique(prop_obs$index)
    
    
    # create log periodograms for replications in ith time seg and jth cov seg
    nseg_time_temp=interval_curr$interval[i,]
    tmin=interval_curr$tmin[i,]
    tmax=interval_curr$tmax[i,]
    
    y <- list()
    yy <- list()
    uu=0
    
    for(l in prop_obs_index)
    {
      uu=uu+1
      
      nfreq <- floor(nseg_time_temp[l] / 2)
      
      y[[uu]] <- log(abs(fft(x_t[tmin[l]:tmax[l], l])) ^ 2 / nseg_time_temp[l])
      yy[[uu]] <- y[[uu]][1:(nfreq + 1)]
      
    }
    
    
    # initial beta
    param=rep(0,nbeta)
    
    
    # BFGS
    postbeta_BFGS <- optim(param,fn, gr, method="BFGS", prop_obs_index, interval_curr$interval[i,], yy, tau_prop[i], nbeta, nbasis, sigmasqalpha)
    
    
    ### sampling & density
    Q=he(postbeta_BFGS$par,prop_obs_index, interval_curr$interval[i,], yy, tau_prop[i], nbeta, nbasis, sigmasqalpha)
    
    if(pre_beta){
      
      ## sampling
      beta_prop=Chol_sampling(Q,nbeta,TT$beta[i,])
      ## Density
      log_beta_prop=Chol_density(beta_prop, Q,nbeta,TT$beta[i,])
      
    }else{
      
      ## sampling
      beta_prop=Chol_sampling(Q,nbeta,postbeta_BFGS$par)
      ## Density
      log_beta_prop=Chol_density(beta_prop,Q,nbeta,postbeta_BFGS$par)
    }
    
    
    
    log_beta_prior_prop=log(dmvnorm(t(beta_prop),matrix(0,nbeta,1),diag(c(sigmasqalpha, tau_prop[i]*matrix(1,nbasis,1))),log = F)) # Prior Density of beta
    
    
    fhat_prop=nu_mat%*%matrix(beta_prop,nbeta,1)
    log_prop_spec_dens=whittle_like(y,prop_obs_index,interval_curr$interval[i,],beta_prop,nbasis) # Loglikelihood  at proposed values
    loglike_prop=log_prop_spec_dens
    
    log_proposal_prop = log_beta_prop
    log_proposal_curr = TT$log_beta[i]
    
    
    ## Prior of beta, tau, g
    lprior_prop=log_beta_prior_prop
    lprior_curr=TT$beta_lprior[i]
    
    
    
    # No boundary cases to correct for with a move step
    
    # log-MH ratio
    alpha[i]=min(1,exp(loglike_prop+lprior_prop-TT$loglike[i]-lprior_curr+log_proposal_curr-log_proposal_prop))
    
    
    rr=runif(1)
    if(alpha[i]>rr)
    {
      
      TT$loglike[i] = loglike_prop
      TT$log_beta[i]=log_beta_prop
      TT$beta_lprior[i]=log_beta_prior_prop
      TT$beta[i,]=beta_prop
      TT$fhat_prop[,i]=fhat_prop
      
      
      accepted[i] = 1
      
      
    }else{
      
      
      alpha[i]=-1
      accepted[i] = 0
      
      
    }
    
    
    
  }
  
  list(alpha=alpha, TT=TT, accepted=accepted)
  
}
