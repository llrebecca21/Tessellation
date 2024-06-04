posterior_eta_br = function(sumPsi, br, Psi, y, eta_r, Sigma){
  -crossprod(sumPsi, br) - sum(y / exp(Psi %*% br)) - sum(c(1,eta_r) * br * br / Sigma)
}