### Whittle likelihood for each tessellation
# y is a list of log periodogram for ith tessellation and lth time series with T[i,l] length, T[i,l] is the ith tessellation for lth time series
# y is y[i][[l]]
# fhat is a list of log power spectrum for ith tessellation and lth time series with T[i,l] length
# beta is the matrix of collection of basis parameters for each tessellation beta[i,j]: ith tessellation, jth basis 

#### whittle_like is for ith tessellation: T is T[i,] and beta[i,]

# source outside files
source("R/lin_basis_func.R")

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
