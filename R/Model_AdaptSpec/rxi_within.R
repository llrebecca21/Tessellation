#get proposal xi

rxi_within = function(tmin, xi_cur,w = 0.3){
  S_cur = length(xi_cur)-1
  # sample m_star wth equal probabilities
  m_star = sample(seq(1,S_cur - 1), 1)
  # valid t's that could be xi_p
  t_seq = seq(xi_cur[m_star] + tmin, xi_cur[m_star+2]-tmin)
  short_check = as.numeric(abs(t_seq - xi_cur[m_star+1]) <= 1)
  # normalize short_check
  short_check = short_check/sum(short_check)
  long_check = rep(1/length(t_seq),length(t_seq))
  #long_move = sample(t_seq,1,long_check)
  ls_check = w*(long_check) + (1-w)*(short_check)
  xi_p = t_seq[sample(1:length(t_seq),1,prob= ls_check)]
  # get xi_prop
  xi_prop = xi_cur
  xi_prop[m_star + 1] = xi_p
  return(list("xi_prop" = xi_prop, "m_star" = m_star))
  
  
}

# xi_prop = rxi_within(tmin = tmin, xi_cur = xi_cur, w = 0.2)$xi_prop
# dxi_within(xi_prop, tmin = tmin, xi_cur = xi_cur, w = 0.2)
