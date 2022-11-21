# functions called inside fn:
#       lin_basis_func
fn = function(param,
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
