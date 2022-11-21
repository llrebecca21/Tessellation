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
