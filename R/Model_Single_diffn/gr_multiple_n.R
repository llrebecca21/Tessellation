gr_multiple_n = function(b, Psi_list, sumsumPsi, perio_list, R, Sigma){
  rs = 0
  for(r in 1:R){
    rs = rs + colSums(Psi_list[[r]] * c(perio_list[[r]] / exp(Psi_list[[r]] %*% b)))
  }
  return(-b/Sigma - sumsumPsi + rs)
}