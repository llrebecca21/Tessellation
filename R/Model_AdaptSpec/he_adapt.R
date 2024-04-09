# Function that calculates the Hessian of the conditional posterior of the beta's
he_adapt = function(Psi, sumPsi, b, Sigma, perio){
  he_bb = -0.5 * crossprod(Psi , Psi * c(perio / exp(Psi %*% b)))
  diag(he_bb) =  -1 / Sigma + diag(he_bb)
  return(he_bb)
}

