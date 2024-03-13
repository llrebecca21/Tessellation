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
source("R/Chol_sampling.R") # performs Cholesky Sampling
source("R/arma_spec.R") # determines the spectral density of the resulting stationary segment
source("R/Model_Single/mcmc_stationary.R") # performs mcmc sampling
source("R/Data_Generation/data_generation.R") # data generating code

# Create a single non-stationary time series
# call function generate_nonstat_abrupt from data_generation.R file

# Determine current number of segments to begin with s^c

# Determine the starting segment break placements \xi_m based on number of segments and tmin condition

#############
# MCMC Moves
# contained in a for loop until convergence

####################
# Within-Model Moves
# Create a function to store this type of move
# Steps:
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








