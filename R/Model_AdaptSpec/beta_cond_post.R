beta_cond_post = function(b, Psi,sumPsi, perio, Sigma){
  -1/2 * (sum(b * b / Sigma) + sum(sumPsi*b) + sum(perio/exp(Psi %*% b)))
}