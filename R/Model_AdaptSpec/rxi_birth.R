rxi_birth = function(Smax, tmin, xi_cur){
  S_cur = length(xi_cur)-1
  len_s = diff(xi_cur)
  # check how many current segments contain at least 2*tmin values
  check_tmin = max(len_s) >= (2*tmin)
  # Find which segment can be chosen
  w = which(len_s >= 2*tmin)
  m_star = w[sample(1:length(w),1)]
  # xi_p
  i = (xi_cur[m_star]+tmin):(xi_cur[m_star+1]-tmin)
  xi_p = i[sample(length(i),1)]
  # insert xi_p into xi_p
  xi_prop = append(xi_cur,xi_p,m_star)
  return(xi_prop)
}







