# Single Stationary Time Series

set.seed(1080)

# Create single time series: AR(1)

# set hyper-parameters

# length of time series
n <- 1000
# burn-in period
burn <- 50
# Create coefficient phi
phi <- 0.5

# Model 1 : AR(1) with phi = 0.5
ts1 <- arima.sim(model = list("ar" = phi), n = n, n.start = burn)






