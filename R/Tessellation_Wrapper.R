# main R script
# Cleaning the main_11-12 file

# clear environment
rm(list = ls())

# Read in the files that are needed
source('R/M_lprior.R')
source('R/time_interval.R')
source('R/range01.R')
Rcpp::sourceCpp('src/distance_partitionC_Lee.cpp')

# Call libraries
library(RcppArmadillo)
library(microbenchmark)

###################################
# Abrupt-abrupt simulate data
#################################
set.seed(1080)

# Define variables:
Ntime=1000
NumX=2
NumXT=NumX+1
NumObs=20

## create covariates matrix of x: dim: 20 x 2
xx=matrix(data = 0, nrow = NumObs, ncol = NumX)

# Question what is True S stand for?
## True S is 2,8,13,18
xx[1,]=c(0.2,0.4)
xx[2,]=c(0.3,0.3)
xx[3,]=c(0.4,0.2)
xx[4,]=c(0.2,0.2)
xx[5,]=c(0.4,0.4)
xx[6,]=c(0.2,0.6)  
xx[7,]=c(0.2,0.8)
xx[8,]=c(0.3,0.7)  
xx[9,]=c(0.4,0.6)  
xx[10,]=c(0.4,0.8)  
xx[11,]=c(0.6,0.2)  
xx[12,]=c(0.6,0.4)  
xx[13,]=c(0.7,0.3)
xx[14,]=c(0.8,0.2)
xx[15,]=c(0.8,0.4)
xx[16,]=c(0.6,0.6)
xx[17,]=c(0.6,0.8)
xx[18,]=c(0.7,0.7)
xx[19,]=c(0.8,0.6)
xx[20,]=c(0.8,0.8)

# Redid code to create column names (got rid of for loop) #
x = matrix(data = 0, nrow = Ntime*NumObs, ncol = NumXT+2)
# Create column names for x1 and x2
xname = paste("x", 1:(NumX), sep = "")
# Create name vector
names = append(c("index", "time", "scaled_time"),xname)
# change matrix x into a dataframe
x = data.frame(x)
# rename the columns of x
colnames(x) = names

# Fill in the columns of x with data #
# replace index with 1:20 with each repeated 1_000 times
x$index=rep(1:NumObs,each=Ntime)
# replace time with 1:1000 repeated 20 times
x$time=rep(1:Ntime,NumObs)
# replace scaled_time with proportion of length completed; repeats each 1_000 index  
x$scaled_time=range01(x$time)

# Repeat each row of xx 1000 times in x
x[,4:(NumXT+2)]=apply(X = xx, MARGIN = 2, FUN = function(i) rep(i,each=Ntime))

## create time series x_t, abrupt change at x2, and at half time points
## Redid ts.sim (unneeded), ts.sim1, ts.sim2 ##
ts.sim1 = matrix(data = 0, nrow = Ntime/2, ncol = NumObs)
ts.sim2 = matrix(data = 0, nrow = Ntime/2, ncol = NumObs)



for(i in 1:NumObs)
{
  if(xx[i,1]<0.5)
  {
    
    if(xx[i,2]<0.5)
    {
      ts.sim1[,i]=arima.sim(list(order=c(1,0,0),ar=0.3),n=Ntime/2)
      ts.sim2[,i]=arima.sim(list(order=c(1,0,0),ar=-0.5),n=Ntime/2)
      
    }else{
      
      ts.sim1[,i]=arima.sim(list(order=c(1,0,0),ar=-0.9),n=Ntime/2)
      ts.sim2[,i]=arima.sim(list(order=c(1,0,0),ar=0.7),n=Ntime/2)
      
    }
    
    
  }else{
    
    if(xx[i,2]<0.5)
    {
      
      ts.sim1[,i]=arima.sim(list(order=c(1,0,0),ar=-0.3),n=Ntime/2)
      ts.sim2[,i]=arima.sim(list(order=c(1,0,0),ar=0.5),n=Ntime/2)
      
    }else{
      
      ts.sim1[,i]=arima.sim(list(order=c(1,0,0),ar=0.9),n=Ntime/2)
      ts.sim2[,i]=arima.sim(list(order=c(1,0,0),ar=-0.7),n=Ntime/2)
      
    }
    
  }
}


x_t=rbind(ts.sim1,ts.sim2)  


# Redid Standardization of x_t #
xmat = cbind(1,1:Ntime)
linfit = solve(crossprod(xmat), crossprod(xmat,x_t))
x_t = x_t - xmat %*% linfit


# Check simulated time series
ts.plot(x_t[,20])


plot(x[,c(4,5)])


###############################################
#  Hyperparameters for tessellation partition
###############################################
# set the seed for the section
set.seed(2333)

## Initialization of tessellation

# Define variables for the tessellation
# number of centers
M=8 
# center locations
S=c(1250,1750,7250,7750,12250,12750,17250,17750) 
# weights
w=c(0.04190068, 0.63559150, 0.32250782) 

# call distance_partitionC function (built in C++)
prt=distance_partitionC_Lee(as.matrix(x[,-c(1,2)]),S,w)

# microbenchmark(distance_partitionC_Lee(as.matrix(x[,-c(1,2)]),S,w), times = 100)

# Call time_interval function (built in R)
interval_curr=time_interval(x,S,prt,NumObs)



















