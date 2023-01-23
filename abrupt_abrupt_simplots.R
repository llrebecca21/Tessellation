# Abrupt-Abrupt Time Series Plots

# First above, second above
ts.plot(x_t[,20])
ts.plot(x_t[,19])
ts.plot(x_t[,18])
ts.plot(x_t[,17])
ts.plot(x_t[,16])

# First Above, Second below
ts.plot(x_t[,15])
ts.plot(x_t[,14])
ts.plot(x_t[,13])
ts.plot(x_t[,12])
ts.plot(x_t[,11])

# First below, second above
par(mfrow = c(2,3), mar = c(4.5,2.5,2,.2))
ts.plot(x_t[,10])
ts.plot(x_t[,9])
ts.plot(x_t[,8])
ts.plot(x_t[,7])
ts.plot(x_t[,6])



# both below: AR1(0.3) and AR1(-0.5) 
# ts.sim1[i,] = arima.sim(list(order = c(1,0,0), ar = 0.3), n = Ntime/2)
# ts.sim2[i,] = arima.sim(list(order = c(1,0,0), ar = -0.5), n = Ntime/2)

par(mfrow = c(2,3), mar = c(4.5,2.5,3,.2), las = 1)


ts.plot(x_t[,1])
ts.plot(x_t[,2])
ts.plot(x_t[,3])
ts.plot(x_t[,4])
ts.plot(x_t[,5])

mtext("AR(1), \U03D5 = 0.3 and \U03D5 = -0.5", outer = TRUE, line = -2)














