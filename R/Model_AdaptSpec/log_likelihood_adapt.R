log_likelihood_adapt = function(b, Psi,sumPsi, perio,n){
  return(-1/2 * (sum(sumPsi*b) + sum(perio/exp(Psi %*% b))) - (n/2)*log(2*pi))
}