# Run a comparison between the different sampling algorithms and plot their mean spectral density on the same plot.
set.seed(105)
par(mfrow = c(2,4))
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
# First run with AR(p) data generating function as the input for both Samplers
# phi = 0.5
R = 8
B = 15
n = 1000
#gendata = generate_adapt(phi = phi, n = n, R = R)
#timeseries = gendata$matrix_timeseries

peaks1 <- runif(R, min = 0.19, max = 0.21)
bandwidths1 <- rep(0.04, R)
peaks2 <- runif(R, min = (pi / 4) - 0.01, max = (pi/4) + 0.01)
bandwidths2 <- rep(0.03, R)

gendata1 = generate_ar2_peak(peaks = peaks1, bandwidths = bandwidths1, n = n)
gendata2 = generate_ar2_peak(peaks = peaks2, bandwidths = bandwidths2, n = n)
timeseries = gendata1$matrix_timeseries + gendata2$matrix_timeseries
ts.plot(timeseries[,1])


# run the Sampler_Single function
#Result_Single = Sampler_Single(timeseries = timeseries, B = B)

# run the Sampler_Wishart function
Result_Wishart = Sampler_Wishart(timeseries = timeseries, B = B, tausquared = 1)

# run the Sampler_eta_br function
Result_eta_br = Sampler_eta_br(timeseries = timeseries, B = B, tausquared = 1)

# run the Sampler_eta_r function
Result_eta_r = Sampler_eta_r(timeseries = timeseries, B = B, tausquared = 1)

# Plot the results below
# Define omega and Psi for plotting purposes
J = floor((n-1) / 2)
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4), mar = c(4,4,3,0.1)+0.1)
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  specdens_eta_r = exp(Psi %*% t(Result_eta_r$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,10), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  # Plot Model Single
  for(h in sample(ncol(specdens_Single), 100, replace = FALSE)){ #light purple
    lines(x = omega, y = log(specdens_Single[,h]), col = rgb(.76, .65, .81, 0.4))
  }
  # Plot Model Wishart
  for(h in sample(ncol(specdens_Wishart), 100, replace = FALSE)){ #dark purple
    lines(x = omega, y = log(specdens_Wishart[,h]), col = rgb(.48, .19, .58, 0.4))
  }
  # Plot Model eta_r
  for(h in sample(ncol(specdens_eta_r), 100, replace = FALSE)){ #light green
    lines(x = omega, y = log(specdens_eta_r[,h]), col = rgb(.65, .85, .62, 0.4))
  }
  # Plot Model eta_br
  for(h in sample(ncol(specdens_eta_br), 100, replace = FALSE)){ #dark green
    lines(x = omega, y = log(specdens_eta_br[,h]), col = rgb(0, .53, .21, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, phi = gendata1$phi[,r]) + arma_spec(omega = omega, phi = gendata2$phi[,r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",
                             #Single
                             rgb(.76, .65, .81),
                             #Wishart
                             rgb(.48, .19, .58),
                             #eta_r
                             rgb(.65, .85, .62),
                             #eta_br
                             rgb(0, .53, .21)),
         lwd = c(2,1,1,1, 1), legend = c("True","Single", "Wishart", "eta_r","eta_br"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
# par(mfrow = c(2,4))
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  specdens_eta_r = exp(Psi %*% t(Result_eta_r$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,10), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Wishart)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = gendata1$phi[,r]) + arma_spec(omega = omega, phi = gendata2$phi[,r])),
        col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  legend("topright", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(3,3,3,3,3),
         legend = c("True", "Single","Wishart","eta_r","eta_br"))
}


###########################
# Run with generate_Krafty
###########################
R = 8
B = 10
n = 1000
set.seed(100)
gendata = generate_Krafty(n = n, R = R)
timeseries = gendata$matrix_timeseries
theta_true = gendata$theta_true

# run the Sampler_Single function
#Result_Single = Sampler_Single(timeseries = timeseries, B = B)

# run the Sampler_Wishart function
Result_Wishart = Sampler_Wishart(timeseries = timeseries, B = B)

# run the Sampler_eta_br function
Result_eta_br = Sampler_eta_br(timeseries = timeseries, B = B, tausquared = 1)

# run the Sampler_eta_r function
Result_eta_r = Sampler_eta_r(timeseries = timeseries, B = B, tausquared = 1)


# Plot the results below
# Define omega and Psi for plotting purposes
J = floor((n-1) / 2)
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

#Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4), mar = c(4,4,3,0.1)+0.1)
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  specdens_eta_r = exp(Psi %*% t(Result_eta_r$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  # Plot Model Single
  for(h in sample(ncol(specdens_Single), 100, replace = FALSE)){ #light purple
    lines(x = omega, y = log(specdens_Single[,h]), col = rgb(.76, .65, .81, 0.4))
  }
  # Plot Model Wishart
  for(h in sample(ncol(specdens_Wishart), 100, replace = FALSE)){ #dark purple
    lines(x = omega, y = log(specdens_Wishart[,h]), col = rgb(.48, .19, .58, 0.4))
  }
  # Plot Model eta_r
  for(h in sample(ncol(specdens_eta_r), 100, replace = FALSE)){ #light green
    lines(x = omega, y = log(specdens_eta_r[,h]), col = rgb(.65, .85, .62, 0.4))
  }
  # Plot Model eta_br
  for(h in sample(ncol(specdens_eta_br), 100, replace = FALSE)){ #dark green
    lines(x = omega, y = log(specdens_eta_br[,h]), col = rgb(0, .53, .21, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",
                             #Single
                             rgb(.76, .65, .81),
                             #Wishart
                             rgb(.48, .19, .58),
                             #eta_r
                             rgb(.65, .85, .62),
                             #eta_br
                             rgb(0, .53, .21)),
         lwd = c(2,1,1,1, 1), legend = c("True","Single", "Wishart", "eta_r","eta_br"))
}


# Plot the Posterior Mean Estimates with True Spectral Density
# par(mfrow = c(2,4))
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  specdens_eta_r = exp(Psi %*% t(Result_eta_r$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-1.5,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Wishart)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  legend("topright", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(3,3,3,3,3),
         legend = c("True", "Single","Wishart","eta_r","eta_br"))
}


