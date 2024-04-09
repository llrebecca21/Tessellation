log_likelihood_adapt = function(b, Psi,sumPsi, perio){
  return(-1/2 * (sum(sumPsi*b) + sum(perio/exp(Psi %*% b))))
}