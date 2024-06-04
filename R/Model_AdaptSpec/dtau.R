dtau = function(tau_p1, tau_p2){
  #return(log(tau_cur/(tau_cur + tau_p1)^2))
  return(-log(2*(tau_p1+tau_p2)^2))
}