rxi_death = function(xi_cur){
  return(xi_cur[-sample(seq(2,length(xi_cur)-1), 1)])
}