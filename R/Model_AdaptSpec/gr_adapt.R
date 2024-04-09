# function that calculates the gradient of the conditional posterior of the beta's
gr_adapt = function(Psi, sumPsi, b, Sigma, perio){
  return(-0.5 * (c(sumPsi) + 2*b/Sigma - colSums(Psi * c(perio / exp(Psi %*% b)))))
}