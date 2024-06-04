he_eta_br = function(Psi, y, br, eta_r, Sigma){
  he1 = -crossprod(Psi , Psi * c(y / exp(Psi %*% br)))
  diag(he1) = diag(he1) - 2 * c(1,eta_r) / Sigma
  return(he1)
}