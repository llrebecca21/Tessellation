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
source("R/Sampling_Algorithms/Sampler_Wishart.R")
source("R/Sampling_Algorithms/Sampler_Single.R")
source("R/Sampling_Algorithms/Sampler_eta_br.R")
# First run with AR(p) data generating function as the input for both Samplers
phi = 0.5
R = 8
B = 10
n = 1000
gendata = generate_adapt(phi = phi, n = n, R = R)
timeseries = gendata$matrix_timeseries

# run the Sampler_Single function
#Result_Single = Sampler_Single(timeseries = timeseries, B = B)

# run the Sampler_Wishart function
Result_Wishart = Sampler_Wishart(timeseries = timeseries, B = B, tausquared = 1)

# run the Sampler_eta_br function
Result_eta_br = Sampler_eta_br(timeseries = timeseries, B = B, tausquared = 1)

#plot(Result_Wishart$Theta[,4], type = "l")
#plot(Result_Single$Theta[,2], type = "l")
#image(apply(Result_Wishart$Lambda_array, c(2, 3),mean))
#round(apply(Result_Wishart$Lambda_array, c(2, 3),mean), 4)
#abline(h = 0)

# Plot the results below
# Define omega and Psi for plotting purposes
J = floor((n-1) / 2)
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4), mar = c(4,4,3,0.1)+0.1)
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  for(h in sample(ncol(specdens_Wishart), 100, replace = FALSE)){ #red
    lines(x = omega, y = log(specdens_Wishart[,h]), col = rgb(.88, .29, .19, 0.4))
  }
  for(h in sample(ncol(specdens_Single), 100, replace = FALSE)){ #yellow
    lines(x = omega, y = log(specdens_Single[,h]), col = rgb(1, .71, 0, 0.4))
  }
  for(h in sample(ncol(specdens_eta_br), 100, replace = FALSE)){ #blue
    lines(x = omega, y = log(specdens_eta_br[,h]), col = rgb(.39,.75,.94, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black", rgb(.88, .29, .19),
                             rgb(.39,.75,.94),
                             rgb(1, .71, 0) ),
         lwd = c(2,1,1,1), legend = c("True", "Wishart", "eta_br" ,"Single"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
# par(mfrow = c(2,4))
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-1.5,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Wishart)), col = rgb(.88, .29, .19, 1))
  lines(x = omega, y = rowMeans(log(specdens_Single)), col = rgb(1, .71, 0, 1))
  lines(x = omega, y = rowMeans(log(specdens_eta_br)), col = rgb(.39,.75,.94, 1))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",rgb(.88, .29, .19, 1),
                               rgb(1, .71, 0, 1),
                               rgb(.39,.75,.94, 1)),
         lwd = c(2,1,1,1),
         legend = c("True", "Wishart", "Single", "eta_br"))
}


###########################
# Run with generate_Krafty
###########################
R = 8
B = 10
n = 1000
gendata = generate_Krafty(n = n, R = R)
timeseries = gendata$matrix_timeseries
theta_true = gendata$theta_true

# run the Sampler_Single function
#Result_Single = Sampler_Single(timeseries = timeseries, B = B)

# run the Sampler_Wishart function
Result_Wishart = Sampler_Wishart(timeseries = timeseries, B = B)

# run the Sampler_eta_br function
Result_eta_br = Sampler_eta_br(timeseries = timeseries, B = B, tausquared = 1)

# Plot the results below
# Define omega and Psi for plotting purposes
J = floor((n-1) / 2)
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4))
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  for(h in sample(ncol(specdens_Wishart), 100, replace = FALSE)){
    lines(x = omega, y = log(specdens_Wishart[,h]), col = rgb(.88, .29, .19, 0.4))
  }
  for(h in sample(ncol(specdens_Single), 100, replace = FALSE)){
    lines(x = omega, y = log(specdens_Single[,h]), col = rgb(1, .71, 0, 0.4))
  }
  for(h in sample(ncol(specdens_eta_br), 100, replace = FALSE)){
    lines(x = omega, y = log(specdens_eta_br[,h]), col = rgb(.39,.75,.94, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",rgb(.88, .29, .19, 1),
                             rgb(1, .71, 0, 1),
                             rgb(.39,.75,.94, 1)),
         lwd = c(2,1,1,1),
         legend = c("True", "Wishart", "Single", "eta_br"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
par(mfrow = c(2,4))
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_eta_br = exp(Psi %*% t(Result_eta_br$bb_beta_array[,,r]))
  Result_Single = Sampler_Single(timeseries = timeseries[,r, drop = FALSE], B = B)
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Wishart)), col = rgb(.88, .29, .19, 1))
  lines(x = omega, y = rowMeans(log(specdens_Single)), col = rgb(1, .71, 0, 1))
  lines(x = omega, y = rowMeans(log(specdens_eta_br)), col = rgb(.39,.75,.94, 1))
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",rgb(.88, .29, .19, 1),
                             rgb(1, .71, 0, 1),
                             rgb(.39,.75,.94, 1)),
         lwd = c(2,1,1,1),
         legend = c("True", "Wishart", "Single", "eta_br"))
}








