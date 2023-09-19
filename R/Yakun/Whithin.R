#' Whithin Function
#'
#' @param x 
#' @param x_t 
#' @param TT 
#' @param NumObs 
#' @param nbasis 
#' @param sigmasqalpha 
#' @param pre_beta 
#'
#' @return
#' @export
#'
#' @examples
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
