# Create function for update of beta and log likelihood in R and then transfer to C++
# Source outside functions
# gr ; fn ; he ; Chol_sampling ; Chol_density ; whittle_like
source("R/fn.R")
source("R/gr.R")
source("R/Chol_sampling.R")
source("R/Chol_density.R")
source("R/whittle_like.R")

update_beta_like = function(x, M, prt, ){
  
  # update beta and loglike
  for(i in 1:M)
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
    
    for(d in prop_obs_index)
    {
      uu=uu+1
      
      nfreq <- floor(nseg_time_temp[d] / 2)
      
      y[[uu]] <- log(abs(fft(x_t[tmin[d]:tmax[d], d])) ^ 2 / nseg_time_temp[d])
      yy[[uu]] <- y[[uu]][1:(nfreq + 1)]
      
    }
    
    
    # initial beta
    param=rep(0,nbeta)
    
    
    # BFGS
    postbeta_BFGS <- optim(param,fn, gr, method="BFGS", prop_obs_index, interval_curr$interval[i,], yy, tau[i], nbeta, nbasis, sigmasqalpha)
    
    
    ### sampling
    Q=he(postbeta_BFGS$par,prop_obs_index, interval_curr$interval[i,], yy, tau[i], nbeta, nbasis, sigmasqalpha)
    beta_prop[i,]=Chol_sampling(Q,nbeta,postbeta_BFGS$par)
    
    
    
    ## Density
    log_beta[i]=Chol_density(beta_prop[i,],Q,nbeta,postbeta_BFGS$par)
    
    
    
    ## prior density
    beta_lprior[i]=dmvnorm(beta_prop[i,],matrix(0,nbeta,1),diag(c(sigmasqalpha, tau[i]*matrix(1,nbasis,1))),log = T) # Prior Density of beta
    
    
    # estimated log of power spectrum and Whittlelikelihood
    fhat_prop[,i]=nu_mat%*%matrix(beta_prop[i,],nbeta,1)
    log_prop_spec_dens=whittle_like(y,prop_obs_index,interval_curr$interval[i,],beta_prop[i,],nbasis) # Loglikelihood  at proposed values
    loglike[i]=log_prop_spec_dens
    
    
  }
  
}
