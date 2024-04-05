# Function to determine if birth or death will occur based
# on properties of the randomly chosen segment S

get_s_fun = function(S, Smax, tmin, timeseries){
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
  # choice for birth or death (version of only allowing dimension to change by +1 or -1)
  new_s = s + sample(c(1,-1),1,prob = c(prob_birth,prob_death))
  return(new_s)
}










