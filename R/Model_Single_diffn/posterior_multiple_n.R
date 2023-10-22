posterior_multiple_n = function(b, Psi_list, sumsumPsi ,perio_list, R, Sigma){
  rs = 0
  for(r in 1:R){
    rs = rs + sum(perio_list[[r]] / exp(Psi_list[[r]] %*% b))
  }
  return(-1/2 * sum(b * b / Sigma) - sum(sumsumPsi * b) - rs)
}