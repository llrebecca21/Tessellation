# AdaptSPEC

# Stationary Time Series
set.seed(101)
# number of time points
n = 256


phi <- c(1.4256, -0.7344, 0.1296)
sim_1 <- arima.sim(model = list("ar" = phi), n = n, n.start = 2000)
ts.plot(sim_1)

# Within Model
# start with small number of segments
# how many segments we currently have
m = 5
# define max number of segments
M = 20


# start with equally spaced xi's that fall between 0 and 256
xi = round(seq(0,256,length.out = 6))

# visualize these points on the time series
abline(v = xi, col = "red")

# take the 4th break and move it between the 3rd and the 5th break.
new_xi = xi
# using -1 and +1 is the way of controlling t-min
# define t-min
tmin = floor(n/M)
# Create variables that will choose at random which xi to move
rand_xi <- sample(m-1,1) 
# moves the randomly chosen xi to a new position
new_xi[rand_xi] = sample(seq(xi[rand_xi - 1] + tmin , xi[rand_xi + 1] - tmin, 1), 1)
# See how the xi was moved from the original location
cbind(xi,new_xi)

# Plot the newly moved xi
ts.plot(sim_1)
#abline(v = xi, col = "red")
abline(v = new_xi, col = "blue", lty = 2)



# Write a function that does the steps above
move_xi <- function(xi, tmin, p){
  # take the xi and make a copy
  new_xi = xi
  # randomly choose a xi index
  m = length(xi) - 1 
  rand_xi = sample(m-1 , 1) + 1
  # p is the chosen threshold for doing sampling method 1 or sampling method 2
  # generate a value that will turn on or off the sampling method
  method <- rbinom(n = 1, size = 1, prob = p)
  if(method == 1){
    # move the randomly chosen xi to a new location conditional on tmin
    new_xi[rand_xi] = sample(seq(xi[rand_xi - 1] + tmin , xi[rand_xi + 1] - tmin, 1), 1)
  } else{
    new_xi[rand_xi] = sample(seq(max(xi[rand_xi] - 1, xi[rand_xi - 1] + tmin),
                                 min(xi[rand_xi] + 1, xi[rand_xi + 1] - tmin)), 1)
  }
  return(new_xi)
}

# p = 1 means we always choose method 1. 
# AdaptSPEC suggests setting p = 0.2 
cbind(xi, move_xi(xi, tmin = tmin, p = 0.2))










