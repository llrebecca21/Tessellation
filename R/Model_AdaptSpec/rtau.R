rtau = function(tau_cur){
  u = runif(1)
  tau_p1 = tau_cur * u / (1-u)
  return(tau_p1)
} 