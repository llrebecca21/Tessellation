rm(list = ls())
setwd("~/Tessellation")

# Load required libraries
library(pracma)
library(trust)
library(mvtnorm)
library(MASS)
library(NPflow)
library(plotly)
library(invgamma)
library(gtools)
library(styler)
library(forecast)
library(RcppArmadillo)

source('R/M_lprior.R')
source('R/tessellation_likelihood.R')
source('R/time_interval.R')
# Get updated file!
#source('R/MHstep_BFGS_11-19.R')
source('R/effective_partition.R')
# Error
#Rcpp::sourceCpp('src/distance_partitionC_Lee.cpp')
# Error
#Rcpp::sourceCpp('src/distance_partition.cpp')


range01 <- function(x){(x-min(x))/(max(x)-min(x))}

##################
# Code Directory
##################

# Abrupt-abrupt simulate data : lines ...
# Slowly-abrupt simulate data : lines ...
# Hyperparameters for tessellation partition : lines ...
# MCMC : lines ...
# Plot of Mean of Posterior : lines ...
# Diagnostics : lines ...
# Sum of Squared Distances : lines ...

###################################
# Abrupt-abrupt simulate data
###################################
set.seed(1080)
# Define variables
Ntime = 1000
NumX = 2
NumXT = NumX + 1
NumObs = 20

## create covariate matrix of x
## xx : (20 x 2) matrix
xx = matrix(0,NumObs,NumX)


## True S (centers) is 2,8,13,18
xx[1,] = c(0.2,0.4) # both below
xx[2,] = c(0.3,0.3) # both below
xx[3,] = c(0.4,0.2) # both below
xx[4,] = c(0.2,0.2) # both below
xx[5,] = c(0.4,0.4) # both below
xx[6,] = c(0.2,0.6) # first below, second above
xx[7,] = c(0.2,0.8) # first below, second above
xx[8,] = c(0.3,0.7) # first below, second above 
xx[9,] = c(0.4,0.6) # first below, second above 
xx[10,] = c(0.4,0.8) # first below, second above  
xx[11,] = c(0.6,0.2) # first above, second below  
xx[12,] = c(0.6,0.4) # first above, second below  
xx[13,] = c(0.7,0.3) # first above, second below
xx[14,] = c(0.8,0.2) # first above, second below
xx[15,] = c(0.8,0.4) # first above, second below
xx[16,] = c(0.6,0.6) # first above, second above
xx[17,] = c(0.6,0.8) # first above, second above
xx[18,] = c(0.7,0.7) # first above, second above
xx[19,] = c(0.8,0.6) # first above, second above
xx[20,] = c(0.8,0.8) # first above, second above

#plot(xx, main = "Plot of Covariances")


# Redo code to create column names # 20_000 x 5
x = matrix(data = 0, nrow = Ntime*NumObs, ncol = NumXT+2)
nrow(x) # 20_000
# Create column names for x1 and x2 
xname = paste("x", 1:(NumX), sep = "")
# Create name vector
names = append(c("index", "time", "scaled_time"),xname)
# change matrix x into a dataframe
x = data.frame(x)
# rename the columns of x
colnames(x) = names

# Fill in the columns of x #
# replace index with 1:20 with each repeated 1_000 times
x$index = rep(1:NumObs, each = Ntime)
# replace time with 1:1000 repeated 20 times
x$time = rep(1:Ntime, NumObs)
# replace scaled_time with proportion of length completed; repeats each 1_000 index  
x$scaled_time = range01(x$time)


# Repeat each row of xx 1000 times 
x[,4:(NumXT+2)] = apply(X = xx, MARGIN = 2, FUN = function(i) rep(i,each = Ntime))


## create time series x_t, abrupt change at x2, and at half time points
## Redo ts.sim1, ts.sim2 ##
ts.sim1 = matrix(data = 0, nrow = NumObs, ncol = Ntime/2)
ts.sim2 = matrix(data = 0, nrow = NumObs, ncol = Ntime/2)

for(i in 1:NumObs)
{
  if(xx[i,1] < 0.5)
  {
    
    if(xx[i,2] < 0.5)
    {
      # Both xx1 and xx2 are below 0.5
      ts.sim1[i,] = arima.sim(list(order = c(1,0,0), ar = 0.3), n = Ntime/2)
      ts.sim2[i,] = arima.sim(list(order = c(1,0,0), ar = -0.5), n = Ntime/2)
      
    }else{
      # xx1 is below 0.5, xx2 is above 0.5
      ts.sim1[i,] = arima.sim(list(order = c(1,0,0), ar = -0.9), n = Ntime/2)
      ts.sim2[i,] = arima.sim(list(order = c(1,0,0), ar = 0.7), n = Ntime/2)
      
    }
    
    
  }else{
    
    if(xx[i,2] < 0.5)
    {
      # xx1 is above 0.5, and xx2 is below 0.5
      ts.sim1[i,] = arima.sim(list(order = c(1,0,0), ar = -0.3), n = Ntime/2)
      ts.sim2[i,] = arima.sim(list(order = c(1,0,0), ar = 0.5), n = Ntime/2)
      
    }else{
      # xx1 is above 0.5, and xx2 is above 0.5
      ts.sim1[i,] = arima.sim(list(order = c(1,0,0), ar = 0.9), n = Ntime/2)
      ts.sim2[i,] = arima.sim(list(order = c(1,0,0), ar = -0.7), n = Ntime/2)
      
    }
    
  }
}

# column bind and transpose the resulting two time series created above
x_t = t(cbind(ts.sim1, ts.sim2))  

# Standardized #
xmat = cbind(1,1:Ntime)
linfit = solve(crossprod(xmat), crossprod(xmat,x_t))
x_t = x_t - xmat %*% linfit

# Plot the 20th time series simulated to check code
# All plots are in abrupt_abrupt_simplots.R
#ts.plot(x_t[,20])
