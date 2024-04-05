# xi_cur = seq(0,1000,by = 100)
# xi_prop = c(xi_cur[1:4],350,xi_cur[5:11])
# tmin = 30
# Smax = floor(xi_cur[length(xi_cur)]/tmin)
get_r_fun = function(Smax, tmin, xi_cur, xi_prop){
  S_cur = length(xi_cur)-1
  S_prop = length(xi_prop)-1
  len_s = diff(xi_cur)
  # check how many current segments contain at least 2*tmin values
  check_tmin = max(len_s) >= (2*tmin)
  # Determine if birth can occur
  can_birth = as.numeric((S_cur < Smax) & check_tmin)
  # Determine if death can occur
  can_death = as.numeric(S_cur > 1)
  # Probability of birth
  prob_birth = can_birth/(can_birth + can_death)
  # Probability of death
  prob_death = 1 - prob_birth
  if(S_prop > S_cur){ # birth
    # r1 = prob_birth
    r1 = prob_birth
    # r2 = prob_segment (based on how many >= 2*tmin)
    # How many segments satisfy >=2*tmin
    num_seg = sum(as.numeric(len_s >= 2*tmin))
    r2 = 1/num_seg
    # r3 = prob_location_in_segment (based on length of chosen segment)
    r3 = 1/(len_s[which(xi_cur != xi_prop[-length(xi_prop)])[1]-1]-2*tmin + 1)
    # r = r1 * r2 * r3
    prob_r = r1*r2*r3
    return(prob_r)
  }else{ # death
    # r1 = prob_death
    r1 = prob_death
    # r2 = prob_xi (based on how many "movable" xi's you have)
    r2 = 1/(length(xi_cur)-2)
    # r = r1 * r2
    prob_r = r1*r2
    return(prob_r)
  }
}