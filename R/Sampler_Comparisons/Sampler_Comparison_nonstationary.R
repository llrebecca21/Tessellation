# Wipe Clean and Set Seed -------------------------------------------------
# Sampler Comparison different lengths of n for each replicate non-stationary time series
rm(list = ls())
# Run a comparison between the different sampling algorithms and plot their mean spectral density on the same plot.
set.seed(100)

# Source Code -------------------------------------------------------------
if(TRUE){
  source("R/Data_Generation/data_generation.R")
  source("R/Model_Single_diffn/posterior_multiple_n.R")
  source("R/Model_Single_diffn/gr_multiple_n.R")
  source("R/Model_Single_diffn/he_multiple_n.R")
  source("R/Model_Single/posterior_multiple.R")
  source("R/Model_Single/gr_multiple.R")
  source("R/Model_Single/he_multiple.R")
  source("R/General_Functions/Chol_sampling.R")
  source("R/General_Functions/arma_spec.R")
  source("R/Data_Generation/data_generation.R")
  source("R/Model_Wishart/gradient_hierarch_Lambda.R")
  source("R/Model_Wishart/posterior_hierarch_Lambda.R")
  source("R/Model_Wishart/he_hierarch_Lambda.R")
  source("R/Model_eta_r/gradient_eta_r.R")
  source("R/Model_eta_r/posterior_eta_r.R")
  source("R/Model_eta_r/he_eta_r.R")
  source("R/Model_eta_br/gradient_eta_br.R")
  source("R/Model_eta_br/posterior_eta_br.R")
  source("R/Model_eta_br/he_eta_br.R")
  source("R/Sampling_Algorithms/Sampler_Wishart.R")
  source("R/Sampling_Algorithms/Sampler_Single.R")
  source("R/Sampling_Algorithms/Sampler_eta_br.R")
  source("R/Sampling_Algorithms/Sampler_eta_r.R")
  source("R/Sampling_Algorithms/Sampler_Single_n.R")
  source("R/Sampling_Algorithms/Sampler_Wishart_n.R")
  source("R/Sampling_Algorithms/Sampler_eta_r_n.R")
  source("R/Sampling_Algorithms/Sampler_eta_br_n.R")
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
  source("R/Sampling_Algorithms/Sampler_adaptspec.R")
  
}


# Set Simulation Parameters -----------------------------------------------
R = 1
B = 10
n = 3000
iter = 1000
tmin = 40
burnin = 10

# Generate Data -----------------------------------------------------------
# First run with AR(p) data generating function as the input for both Samplers

# AR(2) choices
#phi1 = c(0.5, -0.4)
#phi2 = c(0.8, -0.2)

# Choose how many known stationary time series segments we generate data for 
# S_true, either can be prompted in function inputs, or can be randomly generated uniformly from [1,...,Smax]
S_true = 2
# Create S_true # of AR parameters, store in a list
# create a list to store the phi parameters
ar_list = vector(mode = "list", length = S_true)

# sample uniformly [-1,1] for the ar parameters
# for(p in 1:S_true){
#   ar_list[[p]] = c(runif(n = 2, min = -1, max = 1))
# }

# OR 

# Manually pick two AR
phi1 = -0.55
phi2 = 0.55
ar_list[[1]] = phi1
ar_list[[2]] = phi2

# plot the omega
full_omega = seq(0, pi, length.out = 100000)

# Plot the first time series
plot(full_omega, log(arma_spec(full_omega, phi = phi1)), type = "l",
     ylab = "log(spectral density)",
     xlab = "omega",
     main = paste0('Time series 1 phi = ', phi1 ,'\n log(Spectral Density)'))

# Plot the second time series
plot(full_omega, log(arma_spec(full_omega, phi = phi2)), type = "l",
     ylab = "log(spectral density)",
     xlab = "omega",
     main = paste0('Time series 2 phi = ',phi2,'\n log(Spectral Density)'))

# Plot the resulting time series generated
# Initialize a list to store each time series of length "n"
ts_list = vector(mode = "list", length = R)
ar_output = vector(mode = "list", length = R)
s_length = vector(mode = "list", length = R)
for(r in 1:R){
  # n = sample(500:1000, 1)
  # n = avg_n
  ar_output[[r]] = generate_nonstat_abrupt(ar_list = ar_list, n = n, R = R, S_true = S_true)
  ts_list[[r]] = ar_output[[r]]$matrix_timeseries
  s_length[[r]] = ar_output[[r]]$s_length
}

# plot the resulting time series
ts.plot(ts_list[[1]])

# For simplicity for this file, since we only have a single non stationary time series
# set ts_list[[1]] as time series
timeseries = ts_list[[1]]

# Where are the true breaks?
c(0,cumsum(s_length[[1]]))
# 0
# 1509
# 3000

# Run AdaptSPEC function --------------------------------------------------
adaptout = Sampler_adaptspec(timeseries = timeseries, s_length = s_length, iter = iter, tmin = tmin, burnin = 100)

# Retrieve Sampler Output -------------------------------------------------
# If burn-in specified NOT specified in the input:
# Retrieve output from sampler
xi_list = adaptout$xi_list
theta_list = adaptout$theta_list
birth_try = adaptout$birth_try
death_try = adaptout$death_try
birth_yes = adaptout$birth_yes
death_yes = adaptout$death_yes
with_yes = adaptout$with_yes

# If burn-in specified, also retrieve:
xi_list_bi = adaptout$xi_list_bi
theta_list_bi = adaptout$theta_list_bi
birth_try_bi = adaptout$birth_try_bi
death_try_bi = adaptout$death_try_bi
birth_yes_bi = adaptout$birth_yes_bi
death_yes_bi = adaptout$death_yes_bi
with_yes_bi = adaptout$with_yes_bi

# Number of Breakpoints ---------------------------------------------------
# Create the list of the number of breakpoints across all iterations
# and remove the initial iteration from no burn-in since it will appear as 0
# This "initial" value DOES NOT impact the simulations
break_list = sapply(xi_list, FUN = function(x) length(x))[-1]
table(break_list)
# with burn-in removed
break_list_bi = sapply(xi_list_bi, FUN = function(x) length(x))
table(break_list_bi)

# Number of Segments ------------------------------------------------------
# Create the list that stores the number of segments
# and remove the initial iteration from no burn-in since it will appear as -1
# This "initial" value DOES NOT impact the simulations
seg_list = sapply(xi_list, FUN = function(x) length(x)-1)[-1]
table(seg_list)
seg_list_bi = sapply(xi_list_bi, FUN = function(x) length(x)-1)
table(seg_list_bi)

# Posterior Probabilities Segments ----------------------------------------
# with burn-in NOT REMOVED
post_prob_seg = table(sapply(xi_list, FUN = function(x) length(x)-1)[-1])/(iter - 1)
post_prob_seg

# with burn-in REMOVED
post_prob_seg_bi = table(sapply(xi_list_bi, FUN = function(x) length(x)-1))/(iter-burnin)
post_prob_seg_bi

# Plots with burn-in NOT REMOVED ------------------------------------------
xis = unlist(xi_list)
xis = xis[(xis > 0) & (xis < n)]
par(mfrow=c(1,1))
hist(xis, breaks = seq(0,n,10), freq = FALSE,
     main = "Density of the\n Location of Breakpoints",
     xlab = "Location of Breakpoint")
abline(v=c(0,cumsum(s_length[[1]])), col = "red")

# Plot the location of the breakpoints when the number of breakpoints is 3, equivalent to the number of segments is 2
# Create a variable that stores the locations of the breakpoints when 3 breakpoints is chosen
xi_list2 = xi_list[unlist(lapply(xi_list, FUN = length))==3]
par(las = 2)
hist(unlist(xi_list2),
     main = "Location of breakpoints\n for S = 2",
     xlab = "Breakpoint Location")
abline(v=c(0,cumsum(s_length[[1]])), col = "red")
legend("topleft", legend = c("Estimated", "True"), col = c("black","red"),lty = 1)

# Plot the location of the breakpoints when the number of breakpoints is 3, but remove the two end breakpoints
xi2 = unlist(lapply(xi_list, FUN = function(a){if(length(a)==3) return(a[2]) else return(NULL)}))
par(las = 2)
hist(xi2,breaks = seq(0,n,10),
     main = "Location of Breakpoints\n for S = 2 with End Breakpoints Removed",
     xlab = "Breakpoint Location")
abline(v=c(0,cumsum(s_length[[1]])), col = "red")
legend("topleft", legend = c("Estimated", "True"), col = c("black","red"),lty = 1)
# dev.off()

# Plot the location of the breakpoints when the number of breakpoints is 4, equivalent to the number of segments is 3
# Create a variable that stores the locations of the breakpoints when 4 breakpoints is chosen

xi3 = unlist(lapply(xi_list, FUN = function(a){if(length(a)==4) return(a[2:3]) else return(NULL)}))
hist(xi3,breaks = seq(0,n,10),
     main = "Location of Breakpoints\n for S = 3; End BP: Removed",
     xlab = "Breakpoint Location")


# Trace Plots BI NOT Removed ----------------------------------------------

#trace plot of the number of breakpoints before burn-in removed
plot(sapply(xi_list,length)[-1], type = "l",
     main = "Trace Plot: Number of Breakpoints",
     xlab = "Iteration",
     xlim = c(0,1000),
     ylab = "Number of Breakpoints")

# trace plot of the number of segments
plot((sapply(xi_list,length)-1)[-1], type = "l",
     main = "Trace Plot: Number of Segments",
     xlab = "Iteration",
     xlim = c(0,1000),
     ylab = "Number of Segments")


# Plots with burn-in REMOVED ----------------------------------------------
# With burn-in REMOVED
xis_bi = unlist(xi_list_bi)
xis_bi = xis_bi[(xis_bi > 0) & (xis_bi < n)]
hist(xis_bi, breaks = seq(0,n,10), freq = FALSE,
     main = "Density of the\n Location of Breakpoints; BI: Removed",
     xlab = "Location of Breakpoint")
abline(v=c(0,cumsum(s_length[[1]])), col = "red")

# Plot the location of the breakpoints when the number of breakpoints is 3, equivalent to the number of segments is 2
# Create a variable that stores the locations of the breakpoints when 3 breakpoints is chosen
xi_list2_bi = xi_list_bi[unlist(lapply(xi_list_bi, FUN = length))==3]
par(las = 2)
hist(unlist(xi_list2_bi),
     main = "Location of breakpoints\n for S = 2; BI: Removed",
     xlab = "Breakpoint Location")
abline(v=c(0,cumsum(s_length[[1]])), col = "red")
legend("topleft", legend = c("Estimated", "True"), col = c("black","red"),lty = 1)

# Plot the location of the breakpoints when the number of breakpoints is 3, but remove the two end breakpoints
xi2_bi = unlist(lapply(xi_list_bi, FUN = function(a){if(length(a)==3) return(a[2]) else return(NULL)}))
par(las = 2)
hist(xi2_bi,breaks = seq(0,n,10),
     main = "Location of Breakpoints\n for S = 2; BI: Removed; End BP: Removed",
     xlab = "Breakpoint Location")
abline(v=c(0,cumsum(s_length[[1]])), col = "red")
legend("topleft", legend = c("Estimated", "True"), col = c("black","red"),lty = 1)

# Plot the location of the breakpoints when the number of breakpoints is 4, equivalent to the number of segments is 3
# Create a variable that stores the locations of the breakpoints when 4 breakpoints is chosen

xi3_bi = unlist(lapply(xi_list_bi, FUN = function(a){if(length(a)==4) return(a[2:3]) else return(NULL)}))
hist(xi3_bi,breaks = seq(0,n,10),
     main = "Location of Breakpoints\n for S = 3; BI: Removed; End BP: Removed",
     xlab = "Breakpoint Location")

# Trace-Plots BI Removed --------------------------------------------------

# trace plot of the number of breakpoints after burn-in removed
plot(sapply(xi_list_bi,length), type = "l",
     main = "Trace Plot: Number of Breakpoints; BI: Removed",
     xlab = "Iteration",
     xlim = c(0,iter-burnin),
     ylab = "Number of Breakpoints")

# trace plot of the number of segments
plot((sapply(xi_list_bi,length)-1), type = "l",
     main = "Trace Plot: Number of Segments; BI: Removed",
     xlab = "Iteration",
     xlim = c(0,iter-burnin),
     ylab = "Number of Segments")




# Acceptance Ratios -------------------------------------------------------
# Acceptance Ratios with burnin removed

# birth only
birth_ratio_bi = birth_yes_bi/birth_try_bi
# death only
death_ratio_bi = death_yes_bi/death_try_bi
# within only
with_ratio_bi = with_yes_bi/(iter - burnin)

# total acceptance ratio of birth and death move
tot_ratio = (birth_yes + death_yes + with_yes)/(2*(iter-1))

birth_try
death_try

# Plot spectral density
n_o = 100
J = floor(n_o/2)
omega = (2 * pi * (1:J)) / n_o
Psi = cbind(1, outer(X = omega, Y = 1:B, FUN = function(x,y){cos(y * x) * sqrt(2)}))
theta_list[[iter/2]]

j = 520
xi_list[[j]]
par(mfrow = c(1, 2))
plot(omega, Psi%*%theta_list[[j]][1,1:(B+1)], type = "l", ylim = c(-2, 2),
     main = 'Segment 1 estimate')
lines(omega,log(arma_spec(omega = omega,phi = -.5)), col = "blue")
lines(omega,log(arma_spec(omega = omega,phi = .5)), col = "red")
legend('bottomright', col = c("black", "red", "blue"), legend = c("Est", "True 2", 'True 1'), lwd = c(1, 1, 1), bty = 'n')

plot(omega, Psi%*%theta_list[[j]][2,1:(B+1)], type = "l", ylim = c(-2, 2),
     main = 'Segment 2 estimate')
lines(omega,log(arma_spec(omega = omega,phi = -.5)), col = "blue")
lines(omega,log(arma_spec(omega = omega,phi = .5)), col = "red")
legend('bottomright', col = c("black", "red", "blue"), legend = c("Est", "True 2", 'True 1'), lwd = c(1, 1, 1), bty = 'n')


# Create a function that takes a time point 
# loop over iter 
# loop over time
# send in the burnin version
make_heat = function(xi_list, theta_list,omega){
  B = ncol(theta_list[[1]])-2
  # define Psi
  Psi = cbind(1, outer(X = omega, Y = 1:B, FUN = function(x,y){cos(y * x) * sqrt(2)}))
  # 1st dim : number of iterations
  iter = length(xi_list)
  # 2nd dim : length of time series
  time_len = xi_list[[2]][length(xi_list[[2]])]
  # 3rd dim : how many frequencies evaluated at
  len_omega = length(omega)
  # make array to store spectral density estimates at each coordinate
  store_spec = array(data = NA, c(iter, time_len, len_omega))
  # LOOP over iter
  for(i in 1:iter){
    #i = 500
    # LOOP over length of time series
    for(k in 1:time_len){
      #k = 500
      # which xi comes after your current time point
      m = which(xi_list[[i]] >= k)[1]-1
      # grab mth beta vector
      b = theta_list[[i]][m,-ncol(theta_list[[i]])]
      # store the spectral density estimate
      store_spec[i,k,] = Psi %*% b
    }
  }
  quant_spec = apply(store_spec, c(2,3), FUN = function(a){quantile(a,c(0,0.25,0.50,0.75,1))})
  mean_spec = apply(store_spec, c(2,3), mean)
  return(list("quant_spec" = quant_spec, "mean_spec" = mean_spec))
}

# run the function with the burn-in versions

heat_mean = make_heat(xi_list_bi, theta_list_bi, omega)$mean_spec
dim(heat_mean)
# plot the mean of the last time point across all iterations and omegas
# T = 1000
plot(y = heat_mean[n,], x = omega, type = "l")
# T = 1
plot(y = heat_mean[1,], x = omega, type = "l")
# T = 500
plot(y = heat_mean[1500,], x = omega, type = "l")

# blue indicates early time points NOT early iterations
# red indicates later time points
# This plot is good because it shows that the blue dominates the curve corresponding to the true density of the 
# first AR(1) section
# The red dominates the curve that is closest to the AR(1) of the second stationary time series
par(mfrow = c(1,1))
plot(x = c(), y = c(), xlim = range(omega), ylim = range(heat_mean),
     main = "Mean log(Spectral Density Estimate)",
     ylab = "Mean log(Spectral Density)",
     xlab = "omega")
legend("top", legend = c("t = 0", "t = n"), col = c("blue", "red"), lty = c(1,1), lwd = 2)
for(h in 1:n){
  lines(x = omega, y = heat_mean[h,], col = rgb(h/n,0,1-(h/n)))
}

# plotting mean spectral density across the time points
# "time varying spectra"
image(x = 1:n, y = omega, z = heat_mean)
# transposed version with omega on x axis and time on y axis
image(x = omega, y = 1:n, z = t(heat_mean))


# plot the quantiles
heat_quant = make_heat(xi_list_bi, theta_list_bi, omega)$quant_spec
dim(heat_quant)
# min:
image(x = 1:n, y = omega, z = heat_quant[1,,], zlim = range(heat_quant))
# 1st quartile:
image(x = 1:n, y = omega, z = heat_quant[2,,], zlim = range(heat_quant))
# 2nd quartile: (median)
image(x = 1:n, y = omega, z = heat_quant[3,,], zlim = range(heat_quant))
# 3rd quartile:
image(x = 1:n, y = omega, z = heat_quant[4,,], zlim = range(heat_quant))
# max:
image(x = 1:n, y = omega, z = heat_quant[5,,], zlim = range(heat_quant))






