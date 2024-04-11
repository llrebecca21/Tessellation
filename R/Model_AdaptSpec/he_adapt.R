# Function that calculates the Hessian of the conditional posterior of the beta's
he_adapt = function(beta_mstar, Psi, sumPsi, Sigma, perio){
  he_bb = -crossprod(Psi , Psi * c(perio / exp(Psi %*% beta_mstar)))
  diag(he_bb) = diag(he_bb) - 1 / Sigma
  return(he_bb)
}

