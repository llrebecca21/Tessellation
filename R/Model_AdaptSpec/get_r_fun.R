get_r_fun = function(s_cur, s_prop, Smax, tmin){
  # returns the probability for birth or death
  # get current length of the segment S
  # get the indexing that contains the current segment S
  len_S = length(timeseries[a,b])
  # Determine if birth can occur
  can_birth = as.numeric((S < Smax) & (len_S > 2*tmin))
  # Determine if death can occur
  can_death = as.numeric((S > 0) & (len_S > 2*tmin))
  # Probability of birth
  prob_birth = can_birth/(can_birth + can_death)
  # Probability of death
  prob_death = 1 - prob_birth
  if(s_prop > s_cur){
    return(prob_birth)
  }else{
    return(prob_death)
  }
}