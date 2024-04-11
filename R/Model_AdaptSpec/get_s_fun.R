# Function to determine if birth or death will occur based
# on properties of the randomly chosen segment S

get_s_fun = function(Smax, tmin,xi_cur){
  S = length(xi_cur)-1
  # check how many current segments contain at least 2*tmin values
  check_tmin = max(diff(xi_cur)) >= (2*tmin)
  # get current length of the segment S
  len_ts = xi_cur[length(xi_cur)]
  # Determine if birth can occur
  can_birth = as.numeric((S < Smax) & check_tmin)
  # Determine if death can occur
  can_death = as.numeric(S > 1)
  # Probability of birth
  prob_birth = can_birth/(can_birth + can_death)
  # Probability of death
  prob_death = 1 - prob_birth
  # choice for birth or death (version of only allowing dimension to change by +1 or -1)
  new_S = S + sample(c(1,-1),1,prob = c(prob_birth,prob_death))
  return(new_S)
}










