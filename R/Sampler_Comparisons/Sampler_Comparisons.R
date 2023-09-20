# Run a comparison between the different sampling algorithms and plot their mean spectral density on the same plot.
set.seed(100)
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

# First run with AR(p) data generating function as the input for both Samplers
phi = 0.5
R = 8
B = 10
n = 1000
gendata = generate_adapt(phi = phi, n = n, R = R)
timeseries = gendata$matrix_timeseries

# run the Sampler_Single function
Result_Single = Sampler_Single(timeseries = timeseries, B = B)

# run the Sampler_Wishart function
Result_Wishart = Sampler_Wishart(timeseries = timeseries, B = B)



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
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-4,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  for(h in sample(ncol(specdens_Wishart), 100, replace = FALSE)){
    lines(x = omega, y = log(specdens_Wishart[,h]), col = rgb(0, 1, 1, 0.2))
  }
  for(h in sample(ncol(specdens_Single), 100, replace = FALSE)){
    lines(x = omega, y = log(specdens_Single[,h]), col = rgb(1, 0, 1, 0.2))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("bottomleft", col = c("black", "cyan","magenta"), lwd = c(2,1,1), legend = c("True", "Wishart", "Single"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
par(mfrow = c(2,4))
for(r in 1:R){
  specdens_Wishart = exp(Psi %*% t(Result_Wishart$bb_beta_array[,,r]))
  specdens_Single = exp(Psi %*% t(Result_Single$Theta[,-(B+2)]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-4,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
    lines(x = omega, y = rowMeans(log(specdens_Wishart)), col = rgb(0, 1, 1, 1))
    lines(x = omega, y = rowMeans(log(specdens_Single)), col = rgb(1, 0, 1, 1))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("bottomleft", col = c("black", "cyan","magenta"), lwd = c(2,1,1), legend = c("True", "Wishart", "Single"))
}











