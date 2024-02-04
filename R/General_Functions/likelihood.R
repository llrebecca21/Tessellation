likelihood = function(sumPsi, Psi, b, y_bar){
  exp(-1/2 * (crossprod(sumPsi, b) + sum(y_bar / exp(Psi %*% b))))
}