rxi_death = function(xi_cur){
  i = seq(2,length(xi_cur)-1)
  return(xi_cur[-i[sample(length(i), 1)]])
}