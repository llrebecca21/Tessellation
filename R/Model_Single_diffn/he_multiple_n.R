he_multiple_n = function(b, Psi_list, perio_list, Sigma, R){
  rs = 0
  for(r in 1:R){
    rs = rs - crossprod(Psi_list[[r]], Psi_list[[r]] * c(perio_list[[r]] / exp(Psi_list[[r]] %*% b)))
  }
  diag(rs) = diag(rs) - (1 / (Sigma))
  return(rs)
}