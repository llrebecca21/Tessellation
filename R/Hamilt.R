Hamilt=function(beta_temp,index,tmax,tmin,nseg_time_temp,x_t,tau_temp,nbasis,sigmasqalpha,k,ee)
{
  
  
  nbeta <- nbasis + 1
  beta_out = beta_temp
  
  
  
  gr=beta_gradient(beta_temp,index,tmax,tmin,nseg_time_temp,x_t,tau_temp,nbasis,sigmasqalpha)$gradient
  
  
  
  # generate momentum variable
  
  m = matrix(rmvnorm(1,rep(0,nbeta),k*diag(nbeta)),nbeta,1)
  
  #determine leap number and step
  stepsize = runif(1,0,2*ee)
  leap = sample(1:2*ceiling(1/ee),1)
  #leap=5
  
  m_out = m - 0.5*gr*stepsize
  for (i in 1:leap)
  {
    beta_out = beta_out + stepsize*(1/k)*diag(nbeta)%*%m_out
    beta_temp=beta_out
    
    if(i==leap)
    {
      gr=beta_gradient(beta_temp,index,tmax,tmin,nseg_time_temp,x_t,tau_temp,nbasis,sigmasqalpha)$gradient
      cat("gr is ",gr,"\n")
      m_out = m_out - 0.5*gr*stepsize
      
    }else{
      
      gr=beta_gradient(beta_temp,index,tmax,tmin,nseg_time_temp,x_t,tau_temp,nbasis,sigmasqalpha)$gradient
      cat("gr is ",gr,"\n")
      m_out = m_out - 1*gr*stepsize
      
    }
  }
  
  m_out = -m_out
  
  list(beta_out=beta_out,m_out=m_out,m=m)
  
}
