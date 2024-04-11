# AdaptSPEC for Single Non-stationary Time Series
# Uses mcmcstationary.R to estimate each piece wise stationary time series from the non-stationary setup
# Code will be set to R = 1 for now
# This file just runs the simulation 1 time

# libraries:

# Source Code:
set.seed(101)
# Code for updating basis function coefficients for each piece-wise stationary segment
source("R/Model_Single/posterior_multiple.R") # calculates posterior 
source("R/Model_Single/gr_multiple.R") # calculates the gradient
source("R/Model_Single/he_multiple.R") # calculates the hessian
source("R/General_Functions/Chol_sampling.R") # performs Cholesky Sampling
source("R/General_Functions/arma_spec.R") # determines the spectral density of the resulting stationary segment
source("R/Data_Generation/data_generation.R") # data generating code
source("R/Sampling_Algorithms/Sampler_Single.R") # Performs the Sampling algorithm for a stationary time series
source("R/Model_AdaptSpec/get_s_fun.R")
source("R/Model_AdaptSpec/birth_fun.R")
source("R/Model_AdaptSpec/death_fun.R")
source("R/Model_AdaptSpec/within_fun.R")
source("R/Model_AdaptSpec/rxi_death.R")
source("R/Model_AdaptSpec/betapar.R")
source("R/Model_AdaptSpec/get_r_fun.R")
source("R/Model_AdaptSpec/gr_adapt.R")
source("R/Model_AdaptSpec/he_adapt.R")
source("R/Model_AdaptSpec/rxi_birth.R")
source("R/Model_AdaptSpec/dxi_within.R")
source("R/Model_AdaptSpec/rxi_within.R")
source("R/Model_AdaptSpec/dxi.R")
source("R/Model_AdaptSpec/log_likelihood_adapt.R")
source("R/Model_AdaptSpec/beta_cond_post.R")
source("R/Model_AdaptSpec/dtau.R")
source("R/Model_AdaptSpec/rtau.R")
# # Create a single non-stationary time series
# set hyper-parameters
# length of a single time series (entire time-series). Will be prompted in function
n = 1000
# tmin constraint: Will be prompted in function
tmin = 50
# Determine Smax contstraint. Check will occur in function to make sure Smax \ne 0 or 1
Smax = floor(n/tmin) 
# 
R = 1

# Choose how many known stationary time series segments we generate data for 
# S_true, either can be prompted in function inputs, or can be randomly generated uniformly from [1,...,Smax]
S_true = 4
# Create S_true # of AR parameters, store in a list
# sample uniformly [-1,1] for the ar parameters
# create a list to store the phi parameters
ar_list = vector(mode = "list", length = S_true)
for(p in 1:S_true){
  ar_list[[p]] = c(runif(n = 2, min = -1, max = 1))
}

# Plot the resulting time series generated
# Initialize a list to store each time series of length "n"
ts_list = vector(mode = "list", length = R)
ar_output = vector(mode = "list", length = R)
s_length = vector(mode = "list", length = R)
for(r in 1:R){
  # n = sample(500:1000, 1)
  # n = avg_n
  ar_output[[r]] = generate_nonstat_abrupt(ar_list = ar_list, n = n, R = 1, S_true = S_true)
  ts_list[[r]] = ar_output[[r]]$matrix_timeseries
  s_length[[r]] = ar_output[[r]]$s_length
}

# plot the resulting time series
ts.plot(ts_list[[1]])

# For simplicity for this file set ts_list[[1]] as timeseries
timeseries = ts_list[[1]]


# highest little j index value for the frequencies
# J = floor((n-1) / 2)
# Frequency (\omega_j): defined on [0, 2\pi)
# for j = 0,...,n-1
# omega = (2 * pi * (0:(n-1)))/n
# omega = (2 * pi * (0:J)) / n

# Determine the starting segment break placements \xi_m based on number of segments and tmin condition
xi_cur = c(0,cumsum(s_length[[1]]))
# xi_0 = 0 
# xi_1 = 261
# xi_2 = 510
# xi_3 = 760
# xi_4 = 1000

# Specify the intervals for the current time series
# Segment 1: [0, 264]
# Segment 2: [265, 517]
# Segment 3: [518, 748]
# Segment 4: [749, 1000]

B = 10
D = 1 / (4 * pi * (1:B)^2)


# Determine current number of segments to begin with s^c
# To begin set it equal to the truth
S_c = length(xi_cur) - 1
iter = 100
# Initialize storage for the basis coefficients
theta_list = vector(mode = "list", length = iter)
xi_list = vector(mode = "list", length = iter)

# Run the Sampler_Single on each segment to get the initial values for the basis coefficients
# create matrix to store the outputs as matrices
theta_list[[1]] = matrix(NA, nrow = S_c, ncol = B + 2)
for(k in 1:S_c){
  # Run the sampler on the kth segment of the time-series
  seg_cur = timeseries[(xi_cur[k] + 1):(xi_cur[k + 1])]
  # run the sampler and store
  init_out = Sampler_Single(cbind(seg_cur))
  theta_list[[1]][k,] = init_out$Theta[nrow(init_out$Theta),]
}



# Set the current Theta and perio as:
Beta = theta_list[[1]][,-ncol(theta_list[[1]])]
tau = theta_list[[1]][,ncol(theta_list[[1]])]

for(i in 2:iter){
  # Get new S (Proposes a birth or death)
  S_p = get_s_fun(Smax = Smax, tmin = tmin, xi_cur = xi_cur)
  # Compare S_c with S_p
  if(S_p > S_c){ # birth move
    print('DOING BIRTH!')
    output = birth_fun(xi_cur = xi_cur, tmin = tmin, Smax = Smax, Beta = Beta, tau = tau, timeseries = timeseries, sigmasalpha = 100, D = D, B = B, nu_0 = 3)
  }else{
    print('DOING DEATH!')
    output = death_fun(xi_cur = xi_cur, tmin = tmin, Smax = Smax, Beta = Beta, tau = tau, timeseries = timeseries, sigmasalpha = 100, D = D, B = B, nu_0 = 3)
  }
  # Update parameters based on outputs from either the birth or death that was performed
  # update Beta
  Beta = output$Beta
  tau = output$tau
  xi_cur = output$xi
  S_c = length(xi_cur)-1
  
  
  # Within model move
  if(S_c > 1){
    
    print("DOING WITHIN!")
    output = within_fun(xi_cur = xi_cur, tmin = tmin, Smax = Smax, Beta = Beta, tau = tau, timeseries = timeseries, sigmasalpha = 100, D = D, B = B, nu_0 = 3)
    Beta = output$Beta
    tau = output$tau
    xi_cur = output$xi
    
  }
  
  #output to theta_list
  theta_list[[i]] = cbind(Beta,tau)
  #output to xi_list
  xi_list[[i]] = xi_cur

}
###################
# MCMC Moves
# contained in a for loop until convergence

####################
# Within-Model Moves
# within_model_move contains this step
# Steps in the function:

  # 1. given current number of segments, propose a partition point xi_m to be relocated
  
  # 2. The basis functions for the adjacent segments affected by this relocation are updated using mcmcstationary.R 
  
  # Both steps are accepted or rejected jointly in an M-H step
  # Smoothing parameters tau^2 are then updated in a Gibbs Step

#####################
# Between-Model Moves
# Create a function for each type of process

# Steps:

# 1. Determine if move will be birth or death

# 1a. (Birth): If birth move is proposed:
# segment is chosen to split.
# new partition point selected from within this chosen segment
# 2 new smoothing parameters are formed from the current single smoothing parameter
# conditional on these two smoothing parameters, new basis function coefficients are proposed from mcmcstationary.R


# 1b. (Death): If death move is proposed:
# partition point is selected to be removed.
# new smoothing parameter is formed from the previous 2 smoothing parameters
# conditional on the new smoothing parameter, new set of basis function coefficients are proposed.








