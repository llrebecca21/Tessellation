# function that calculates the gradient of the conditional posterior of the beta's
gr_adapt = function(b, Psi, sumPsi, Sigma, perio){
  return(c(sumPsi) + b/Sigma - colSums(Psi * c(perio / exp(Psi %*% b))))
}