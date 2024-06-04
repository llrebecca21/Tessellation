beta_cond_post = function(b, Psi,sumPsi, perio, Sigma,n){
  # This is in log scale already
  1/2 * sum(b * b / Sigma) + sum(sumPsi*b) + sum(perio/exp(Psi %*% b)) + (n/2)*log(2*pi)
}