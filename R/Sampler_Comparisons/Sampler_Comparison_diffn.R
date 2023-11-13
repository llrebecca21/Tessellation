# Sampler Comparison different lengths of n for each replicate time series

# Run a comparison between the different sampling algorithms and plot their mean spectral density on the same plot.
set.seed(105)
par(mfrow = c(2,4))
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
# First run with AR(p) data generating function as the input for both Samplers
phi = c(0.5, -0.4)
#phi = 0.5
R = 8
B = 10
avg_n = 750


# Initialize a list to store each time series of length "n"
ts_list = vector(mode = "list", length = R)
for(r in 1:R){
  # n = rpois(n = 1, lambda = avg_n)
  n = sample(500:1000, 1)
  #n = avg_n
  ts_list[[r]] = generate_adapt(phi = phi, n = n, R = 1)$matrix_timeseries
}

list_n = sapply(ts_list, nrow)
#timeseries = matrix(unlist(ts_list), ncol = R)
list_n


# run the Sampler_eta_br_n function
Result_eta_br_n = Sampler_eta_br_n(ts_list = ts_list, B = B, tausquared = 1)

# run the Sampler_eta_r_n function
Result_eta_r_n = Sampler_eta_r_n(ts_list = ts_list, B = B, tausquared = 1)

# run the Sampler_single_n function
Result_Single_n = Sampler_Single_n(ts_list = ts_list, B = B)

# run the Sampler_Wishart_n function
# V is the precision matrix - hyperparameter for Lambda^{-1}
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(1,1 / (4 * pi * (1:B)^2))))


# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor((n-1) / 2)
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4), mar = c(4,4,3,0.1)+0.1)
for(r in 1:R){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  # Plot Model Single n
  for(h in sample(ncol(specdens_Single_n), 100, replace = FALSE)){ #light purple
    lines(x = omega, y = log(specdens_Single_n[,h]), col = rgb(.76, .65, .81, 0.4))
  }
  # Plot Independent Model n
  for(h in sample(ncol(specdens_Independent_n), 100, replace = FALSE)){ #hotpink
    lines(x = omega, y = log(specdens_Independent_n[,h]), col = rgb(1, .412, .71, 0.4))
  }
  # Plot Model Wishart n
  for(h in sample(ncol(specdens_Wishart_n), 100, replace = FALSE)){ #dark purple
    lines(x = omega, y = log(specdens_Wishart_n[,h]), col = rgb(.48, .19, .58, 0.4))
  }
  # Plot Model eta_r
  for(h in sample(ncol(specdens_eta_r_n), 100, replace = FALSE)){ #light green
    lines(x = omega, y = log(specdens_eta_r_n[,h]), col = rgb(.65, .85, .62, 0.4))
  }
  # Plot Model eta_br_n
  for(h in sample(ncol(specdens_eta_br_n), 100, replace = FALSE)){ #dark green
    lines(x = omega, y = log(specdens_eta_br_n[,h]), col = rgb(0, .53, .21, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",
                             #Single light purple
                             rgb(.76, .65, .81),
                             #Independent hot pink
                             rgb(1, .412, .71, 0.4),
                             #Wishart dark purple
                             rgb(.48, .19, .58),
                             #eta_r light green
                             rgb(.65, .85, .62),
                             #eta_br dark green
                             rgb(0, .53, .21)),
         lwd = c(2,1,1,1, 1), legend = c("True", "Single", "Independent","Wishart", "eta_r","eta_br"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
# par(mfrow = c(2,4))
for(r in 1:R){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-1.5,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  legend("topright", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(1, .412, .71),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(3,3,3,3,3),
         legend = c("True", "Single", "Independent","Wishart","eta_r","eta_br"))
}


###########################
# Run with generate_Krafty
###########################
R = 8
B = 10
set.seed(100)

# Initialize a list to store each time series of length "n"
ts_list = vector(mode = "list", length = R)
theta_true = rep(NA, length = R)
for(r in 1:R){
  # n = rpois(n = 1, lambda = avg_n)
  n = sample(500:1000, 1)
  #n = avg_n
  gendata = generate_Krafty(n = n, R = 1)
  ts_list[[r]] = gendata$matrix_timeseries
  theta_true[r] = gendata$theta_true
}

list_n = sapply(ts_list, nrow)
list_n


#gendata = generate_Krafty(n = n, R = R)
#timeseries = gendata$matrix_timeseries
#theta_true = gendata$theta_true

# run the Sampler_Single_n function
Result_Single_n = Sampler_Single_n(ts_list = ts_list, B = B)

# run the Sampler_Wishart function
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B)

# run the Sampler_eta_br function
Result_eta_br_n = Sampler_eta_br_n(ts_list = ts_list, B = B, tausquared = 1)

# run the Sampler_eta_r function
Result_eta_r_n = Sampler_eta_r_n(ts_list = ts_list, B = B, tausquared = 1)


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
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  # Plot Model Single n
  for(h in sample(ncol(specdens_Single_n), 100, replace = FALSE)){ #light purple
    lines(x = omega, y = log(specdens_Single_n[,h]), col = rgb(.76, .65, .81, 0.4))
  }
  # Plot Independent Model n
  for(h in sample(ncol(specdens_Independent_n), 100, replace = FALSE)){ #hotpink
    lines(x = omega, y = log(specdens_Independent_n[,h]), col = rgb(1, .412, .71, 0.4))
  }
  # Plot Model Wishart n
  for(h in sample(ncol(specdens_Wishart_n), 100, replace = FALSE)){ #dark purple
    lines(x = omega, y = log(specdens_Wishart_n[,h]), col = rgb(.48, .19, .58, 0.4))
  }
  # Plot Model eta_r
  for(h in sample(ncol(specdens_eta_r_n), 100, replace = FALSE)){ #light green
    lines(x = omega, y = log(specdens_eta_r_n[,h]), col = rgb(.65, .85, .62, 0.4))
  }
  # Plot Model eta_br_n
  for(h in sample(ncol(specdens_eta_br_n), 100, replace = FALSE)){ #dark green
    lines(x = omega, y = log(specdens_eta_br_n[,h]), col = rgb(0, .53, .21, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",
                             #Single light purple
                             rgb(.76, .65, .81),
                             #Independent hot pink
                             rgb(1, .412, .71, 0.4),
                             #Wishart dark purple
                             rgb(.48, .19, .58),
                             #eta_r light green
                             rgb(.65, .85, .62),
                             #eta_br dark green
                             rgb(0, .53, .21)),
         lwd = c(2,1,1,1, 1), legend = c("True","Single", "Independent" ,"Wishart", "eta_r","eta_br"))
}


# Plot the Posterior Mean Estimates with True Spectral Density
# par(mfrow = c(2,4))
for(r in 1:R){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-1.5,2), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, theta = theta_true[r])), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  legend("topright", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(1, .412, .71),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(3,3,3,3,3),
         legend = c("True", "Single", "Independent", "Wishart","eta_r","eta_br"))
}




###########################
# Generate an AR(2) Mix
###########################
set.seed(105)
R = 8
B = 25
#avg_n = 750

peaks1 <- runif(R, min = 0.19, max = 0.21)
bandwidths1 <- rep(0.1, R)
peaks2 <- runif(R, min = (pi / 4) - 0.2, max = (pi/4) + 0.2)
bandwidths2 <- rep(0.1, R)

# Initialize a list to store each time series of length "n"
gendata1 = vector(mode = "list", length = R)
gendata2 = vector(mode= "list", length = R)
ts_list = vector(mode = "list", length = R)
for(r in 1:R){
  # n = rpois(n = 1, lambda = avg_n)
  n = sample(500:1000, 1)
  #n = avg_n
  gendata1[[r]] = generate_ar2_peak(peaks = peaks1[r], bandwidths = bandwidths1[r], n = n)
  gendata2[[r]] = generate_ar2_peak(peaks = peaks2[r], bandwidths = bandwidths2[r], n = n)
  ts_list[[r]] = gendata1[[r]]$matrix_timeseries + gendata2[[r]]$matrix_timeseries
  #ts_list[[r]] = generate_ar2_peak(peaks = peaks1, bandwidths = bandwidths1, n = n)$matrix_timeseries
}

list_n = sapply(ts_list, nrow)
#timeseries = matrix(unlist(ts_list), ncol = R)
list_n

ts_list
# run the Sampler_eta_br_n function
Result_eta_br_n = Sampler_eta_br_n(ts_list = ts_list, B = B, tausquared = 1)

# run the Sampler_eta_r_n function
Result_eta_r_n = Sampler_eta_r_n(ts_list = ts_list, B = B, tausquared = 1)

# run the Sampler_single_n function
Result_Single_n = Sampler_Single_n(ts_list = ts_list, B = B)

# run the Sampler_Wishart_n function
# V is the precision matrix - hyperparameter for Lambda^{-1}
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(1,1 / (4 * pi * (1:B)^2))))


# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor((n-1) / 2)
omega = (2 * pi * (0:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,4), mar = c(4,4,3,0.1)+0.1)
for(r in 1:R){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,10), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  # Plot Model Single n
  for(h in sample(ncol(specdens_Single_n), 100, replace = FALSE)){ #light purple
    lines(x = omega, y = log(specdens_Single_n[,h]), col = rgb(.76, .65, .81, 0.4))
  }
  # Plot Independent Model n
  for(h in sample(ncol(specdens_Independent_n), 100, replace = FALSE)){ #hotpink
    lines(x = omega, y = log(specdens_Independent_n[,h]), col = rgb(1, .412, .71, 0.4))
  }
  # Plot Model Wishart n
  for(h in sample(ncol(specdens_Wishart_n), 100, replace = FALSE)){ #dark purple
    lines(x = omega, y = log(specdens_Wishart_n[,h]), col = rgb(.48, .19, .58, 0.4))
  }
  # Plot Model eta_r
  for(h in sample(ncol(specdens_eta_r_n), 100, replace = FALSE)){ #light green
    lines(x = omega, y = log(specdens_eta_r_n[,h]), col = rgb(.65, .85, .62, 0.4))
  }
  # Plot Model eta_br_n
  for(h in sample(ncol(specdens_eta_br_n), 100, replace = FALSE)){ #dark green
    lines(x = omega, y = log(specdens_eta_br_n[,h]), col = rgb(0, .53, .21, 0.4))
  }
  lines(x = omega, y = log(arma_spec(omega = omega, phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topright", col = c("black",
                             #Single light purple
                             rgb(.76, .65, .81),
                             #Independent hot pink
                             rgb(1, .412, .71, 0.4),
                             #Wishart dark purple
                             rgb(.48, .19, .58),
                             #eta_r light green
                             rgb(.65, .85, .62),
                             #eta_br dark green
                             rgb(0, .53, .21)),
         lwd = c(2,1,1,1, 1), legend = c("True", "Single", "Independent","Wishart", "eta_r","eta_br"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
par(mfrow = c(2,4))
for(r in 1:R){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,10), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  legend("topright", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(1, .412, .71),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(3,3,3,3,3),
         legend = c("True", "Single", "Independent","Wishart","eta_r","eta_br"))
}










