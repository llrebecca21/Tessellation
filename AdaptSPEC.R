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














