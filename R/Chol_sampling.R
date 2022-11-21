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
