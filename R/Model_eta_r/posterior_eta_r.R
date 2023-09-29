posterior_eta_r = function(sumPsi, br, Psi, y, eta_r, Sigma){
  -crossprod(sumPsi, br) - sum(y / exp(Psi %*% br)) - sum(eta_r * br * br / Sigma)
}