# AdaptSPEC for Single Non-stationary Time Series
# Uses mcmcstationary.R to estimate each piece wise stationary time series from the non-stationary setup
# Code will be set to R = 1 for now
# This file just runs the simulation 1 time

# libraries:

# Source Code:
set.seed(100)
# Code for updating basis function coefficients for each piece-wise stationary segment
source("R/Model_Single/posterior_multiple.R") # calculates posterior 
source("R/Model_Single/gr_multiple.R") # calculates the gradient
source("R/Model_Single/he_multiple.R") # calculates the hessian
source("R/General_Functions/Chol_sampling.R") # performs Cholesky Sampling
source("R/General_Functions/arma_spec.R") # determines the spectral density of the resulting stationary segment
source("R/Model_Single/mcmc_stationary.R") # performs mcmc sampling
source("R/Data_Generation/data_generation.R") # data generating code

# # Create a single non-stationary time series
# set hyper-parameters
# length of a single time series (entire time-series). Will be prompted in function
n = 1000
# tmin constraint: Will be prompted in function
tmin = 50
# Determine Smax contstraint. Check will occur in function to make sure Smax \ne 0 or 1
Smax = floor(n/tmin) 

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

# highest little j index value for the frequencies
# J = floor((n-1) / 2)
# Frequency (\omega_j): defined on [0, 2\pi)
# for j = 0,...,n-1
# omega = (2 * pi * (0:(n-1)))/n
# omega = (2 * pi * (0:J)) / n

# Determine the starting segment break placements \xi_m based on number of segments and tmin condition
xi_c = c(0,cumsum(s_length[[1]]))

# Determine current number of segments to begin with s^c
# To begin set it equal to the truth
S_c = length(xi_c) - 1

# Run the mcmc_stationary sampler on each segment to get the initial values for the basis coefficients



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








