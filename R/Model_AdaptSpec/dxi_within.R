dxi_within = function(xi_prop, tmin, xi_cur,w = 0.2){
  S_cur = length(xi_cur)-1
  # get m_star that was sampled in the rxi_within function
  m_star = which(xi_cur != xi_prop)-1
  if(length(m_star) == 0){
    
    total_prob = 0
    for(m_star in 1:(S_cur - 1)){
      
      t_seq = seq(xi_cur[m_star] + tmin, xi_cur[m_star+2]-tmin)
      short_check = as.numeric(abs(t_seq - xi_cur[m_star+1]) <= 1)
      # normalize short_check
      short_check = short_check/sum(short_check)
      long_check = rep(1/length(t_seq),length(t_seq))
      #long_move = sample(t_seq,1,long_check)
      ls_check = w*(long_check) + (1-w)*(short_check)
      
      total_prob = total_prob + 1/(S_cur - 1) * ls_check[t_seq == xi_prop[m_star+1]]
      
    }
    return(total_prob)
    
  }
  # valid t's that could be xi_p
  t_seq = seq(xi_cur[m_star] + tmin, xi_cur[m_star+2]-tmin)
  short_check = as.numeric(abs(t_seq - xi_cur[m_star+1]) <= 1)
  # normalize short_check
  short_check = short_check/sum(short_check)
  long_check = rep(1/length(t_seq),length(t_seq))
  #long_move = sample(t_seq,1,long_check)
  ls_check = w*(long_check) + (1-w)*(short_check)
  return(1/(S_cur - 1) * ls_check[t_seq == xi_prop[m_star+1]])
  # 0.03111111
  
  
}

