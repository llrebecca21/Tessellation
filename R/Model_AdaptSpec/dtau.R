dtau = function(tau_p1, tau_cur){
  return(log(tau_cur/(tau_cur + tau_p1)^2))
}