beta_gradient=function(beta_temp,index,tmax,tmin,nseg_time_temp,x_t,tau_temp,nbasis,sigmasqalpha)
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
  
  
  nbeta <- nbasis + 1
  param=beta_temp
  
  
  ## create gradient of log posterior beta probability
  gradient=whittle_derivs2(param,index,nseg_time_temp,yy,tau_temp,nbeta,nbasis,sigmasqalpha)$gradient
  
  list(gradient=gradient,y=y)
  
  
}