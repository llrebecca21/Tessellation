# Within Model Move
within_model_move = function(){
  # Steps:
  # 1. given current number of segments, propose a partition point xi_m to be relocated
  xi_ind = sample(x = seq(2,S_c),size =1)
  
  # 1i. Select a partition point uniformly from the available S_c-1 partition points
  # Determine range new xi_m can be placed
  prop_int_low = xi_c[xi_ind - 1] + tmin
  prop_int_high = xi_c[xi_ind + 1] - tmin
  # Select a proposed value in this range
  # Choose a new spot (small or large step)
  smlr = rbinom(1)
  #if()
  
  xi_p = sample(x = seq(prop_int_low, prop_int_high), size = 1)
  xi_prop = xi_c
  xi_prop[xi_ind] = xi_p
  
  # 2. The basis functions for the adjacent segments affected by this relocation are updated using mcmcstationary.R 
  
  # Both steps are accepted or rejected jointly in an M-H step
  # Smoothing parameters tau^2 are then updated in a Gibbs Step
  
  
}







