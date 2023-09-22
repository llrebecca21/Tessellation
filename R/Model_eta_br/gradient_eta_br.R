gradient_eta_br = function(br, sumPsi, Psi, y, Sigma, eta_r){
  -sumPsi + colSums(Psi * c(y / exp(Psi %*% br))) - 2 * br * eta_r / Sigma
}