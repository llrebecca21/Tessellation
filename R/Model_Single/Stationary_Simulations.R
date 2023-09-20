# Multiple Single Stationary Time Series from the same underlying spectra
library(mvtnorm)
library(progress)
library(fda)
library(bayesplot)
library(ggplot2)
library(mcmcse)
set.seed(100)
source("R/Model_Single/posterior_multiple.R")
source("R/Model_Single/gr_multiple.R")
source("R/Model_Single/he_multiple.R")
source("R/Chol_sampling.R")
source("R/arma_spec.R")
source("R/Model_Single/mcmc_stationary.R")
source("R/Data_Generation/data_generation.R")

# Set outer parameters for simulations
n = 1000
iter = 1000
R = 8
B = 10
phi = 0.5
# this must match what is used in mcmc_stationary.R
burnin = 10

# create vector of frequencies to check
J = floor((n-1) / 2)
test_omega = ((2 * pi * (0:J)) / n)[-1]
omega = (2 * pi * (0:J)) / n
# length(test_omega)
# 249
# True spectral density value at a given omega
true_spec = arma_spec(omega = test_omega, phi = phi, theta = 0)
# Define Psi
test_Psi = outer(X = test_omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
test_Psi[,1] = 1

# Obtaining "true" Betas
full_omega = seq(0, pi, length.out = 100000)
full_spec = arma_spec(omega = full_omega, phi = phi, theta = 0)
full_Psi =  outer(X = full_omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
full_Psi[,1] = 1
true_Beta = coef(lm(log(full_spec) ~ full_Psi + 0))

# How many simulations:
Sim = 1000

# Create object to store coverage checks
Coverage = matrix(NA, nrow = Sim, ncol = length(test_omega))

# Create object to store simultaneous beta checks
Beta_distance = rep(NA, Sim)

# Create object to store bias calculations of beta and of spectral density
# Object to store beta bias:
beta_bias = matrix(NA, nrow = Sim, ncol = B + 1)
# Object to store spectral density bias:
spec_bias = matrix(NA, nrow = Sim, ncol = length(test_omega))
# Store RMSE in a matrix
rmse = matrix(NA, nrow = Sim, ncol = B + 1)

# Create array object to store all of the post_pred samples
post_pred = array(data = NA, dim = c(length(test_omega), iter - burnin, Sim))

# Create vector object to store simultaneous coverage checks for spectral density
Simult_truespec = rep(NA, Sim)

pb = progress_bar$new(total = Sim)
t1 = Sys.time()
# construct for loop for simuation:
for(i in 1:Sim){
  # i = 1
  
  pb$tick()
  
  funccall = mcmc_stationary(n = n, iter = iter, phi = phi, R = R, B = B)
  # extract
  Theta = funccall$Theta
  # dim(Theta)
  av_perio = funccall$av_perio
  
  #
  
  test_specdens = exp(test_Psi %*% t(Theta[ ,-(B+2)]))
  
  # Obtain Posterior Predictive (creates samples from posterior predictive distribution)
  if(TRUE){  
    # As R increases the distribution of post_pred becomes more bell-shaped
    # As n increases it gets closer to true asympototic Gamma.
    # hist(post_pred[1,])
    post_pred[,,i] = matrix(rgamma(nrow(test_specdens) * ncol(test_specdens), R, R / c(test_specdens)),
                            nrow(test_specdens),
                            ncol(test_specdens))
  }

  
  # Simultaneous Coverage for Spectral Density
  if(TRUE){
    Simult_specdens = fbplot(fit = log(test_specdens), x = test_omega, xlim = range(test_omega), ylim = range(log(test_specdens)), prob = c(0.95, 0.5, .05),
                             xlab = "omega",
                             factor = 1.5,
                             barcol = "dodgerblue" ,
                             col = c("purple", "orange", "violetred3"),
                             fullout = TRUE,
                             plot = FALSE)
    deepcurve = log(test_specdens)[, Simult_specdens$depth >= quantile(Simult_specdens$depth, .05)]
    uppercurve = apply(deepcurve, 1, FUN = max)
    lowercurve = apply(deepcurve, 1, FUN = min)
    Simult_truespec[i] = as.numeric(all((log(true_spec) <= uppercurve) & (log(true_spec) >= lowercurve)))
    
  }
  
  
  # Point-wise Spectral Density Coverage
  if(TRUE){
    # Obtain the quantiles from specdens
    crit_bounds = apply(test_specdens, 1, FUN = function(x){quantile(x,c(0.025, 0.975))})
    
    # Check if true_spec is in crit_bounds
    coverage_check = ifelse(true_spec <= crit_bounds[2,] & true_spec >= crit_bounds[1,], yes = 1, no = 0)
    
    # For point-wise for spectral density
    Coverage[i,] = coverage_check
  }

  
  # Simultaneous coverage check for Beta and bias check for Beta
  Beta = Theta[, -c(B+2)]
  
  # Calculate Bias (Beta and Spectral Density): estimate - true
  if(TRUE){
    # Check Bias of Beta's
    beta_bias[i,] = colMeans(Beta) - true_Beta
    
    # Check Bias of Spectral Density
    spec_bias[i,] = rowMeans(test_specdens) - true_spec
    
  }

  # Simultaneous Beta Coverage Check
  if(TRUE){
    # Center Beta
    cent_Beta = Beta - tcrossprod(rep(1,nrow(Beta)), colMeans(Beta))
    # Estimate covariance matrix of Beta (sample covariance matrix)
    Beta_cov = crossprod(cent_Beta) / nrow(cent_Beta)
    # Calculate the weighted distance of true beta and sample beta
    x = (true_Beta - colMeans(Beta))
    Beta_distance[i] = c(t(x) %*% solve(Beta_cov, x))
  }

  
  # RMSE : how far away on average the estimator is from the true value
  if(TRUE){
    rmse[i,] = sqrt(colMeans((Beta - tcrossprod(rep(1,iter - burnin),true_Beta))^2))
  }


}
Sys.time() - t1


# Point-wise Coverage estimate
colMeans(Coverage)
plot.ts(colMeans(Coverage))

# Simultaneous Beta Coverage
print(Beta_distance)
hist(Beta_distance)

# Simultaneous Beta distance critical value
# qchisq(p = .95, df = B + 1)

mean(Beta_distance <= qchisq(p = .95, df = B + 1))
# 0.847

# Beta Bias

par(mfrow = c(3,4))
for(i in 1:ncol(beta_bias)){
  hist(beta_bias[,i],
       main = "Beta Bias",
       xlab = "bias")
}

# Bias
beta_bias_mean = round(colMeans(beta_bias), 4)
beta_bias_sd = apply(beta_bias, 2, sd)
beta_bias_se = beta_bias_sd / sqrt(Sim)
# plot the beta_bias with 2 SE (to visualize)
#pdf(file = "~/Documents/Tessellation/R/Simulation_Results/Beta_bias.pdf")
par(mfrow = c(3,4))
for(i in 1:ncol(beta_bias)){
  hist(beta_bias[,i])
  polygon(x = c(beta_bias_mean[i] - 2 * beta_bias_se[i], 
                beta_bias_mean[i] + 2 * beta_bias_se[i],
                beta_bias_mean[i] + 2 * beta_bias_se[i],
                beta_bias_mean[i] - 2 * beta_bias_se[i]),
          y = c(0, 0, 50, 50), col = "lightblue")
}
#dev.off()
# Do boolean to show which are within 2SE of the mean bias estimate
ifelse((0 >= beta_bias_mean - 2 * beta_bias_se) & (0 <= beta_bias_mean + 2 * beta_bias_se), yes = 1, no = 0)
# Seed 100/R1:
# 0 0 0 0 0 1 0 1 1 1 1
# Seed 31/R5 :
# 1 0 0 0 0 1 1 1 1 1 1
# Seed 100/R5 :
# 1 0 0 0 0 1 1 1 1 1 1
# Seed 31/n1000/R5 :
# 1 0 0 0 1 1 0 1 1 1 1
# Seed 100/n1000/R5 :
# 1 0 0 0 1 1 1 1 1 1 1

############################
# Spectral Bias
############################
# par(mfrow = c(3,4))
# for(i in 1:ncol(spec_bias)){
#   hist(spec_bias[,i])
# }
# Bias
spec_bias_mean = round(colMeans(spec_bias), 4)
spec_bias_sd = apply(spec_bias, 2, sd)
spec_bias_se = spec_bias_sd / sqrt(Sim)

# Plot the mean bias Estimate and the 2 SE 
# par(mfrow = c(1,1))
# plot( x = test_omega, y = spec_bias_mean, type = "l")
# lines(x = test_omega, y = spec_bias_mean + 2 * spec_bias_se, col = "red")
# lines(x = test_omega, y = spec_bias_mean - 2 * spec_bias_se, col = "red")
# abline(h = 0, col = "purple")


# Do boolean to show which are within 2SE of the mean bias estimate
specbias_2se = ifelse((0 >= spec_bias_mean - 2 * spec_bias_se) & (0 <= spec_bias_mean + 2 * spec_bias_se), yes = 1, no = 0)
sum(specbias_2se)
# 64 in
mean(specbias_2se)
# 25% within

# Point-wise (across omegas) Spectral Density Coverage Estimate
par(mfrow = c(1,1))
pdf(file = "~/Documents/Tessellation/R/Simulation_Results/Pointwise_SpecDens_Coverage.pdf")
plot(x = test_omega, y = colMeans(Coverage), type = "l")
dev.off()

#################
# Sim Results
#################
# Inputs: 
# Run Time : 1.573408 mins
# n = 500
# iter = 1000
# R = 1
# B = 10
# Sim = 100
# Simultaneous Coverage (M Distance) :
# 0.82
# Simultaneous Coverage (Fisher's) ;
# 0.84
# Pointwise Coverage :
# 0.70 0.75 0.96 0.91 0.89 0.87 0.88 0.90 0.92 0.93
# Bias Beta :
# -0.012733989 -0.045276903 -0.035237148 -0.011256257 -0.006624983 -0.004984692 -0.003128275 -0.001289735  0.004165166
# -0.002043330 -0.001567824
# Bias Beta Standard Deviation
# 0.06285571 0.06655048 0.05153184 0.04528638 0.03457813 0.03323869 0.03003279 0.02089856 0.01745743 0.01336114 0.01103776
# Bias Beta Standard Error :
# 0.0029217744 0.0029711704 0.0025770661 0.0019541863 0.0016083050 0.0013184106 0.0012095824 0.0008989013
# 0.0008531989 0.0006248765 0.0005571732
# Is 0 within 2 standard errors of the beta bias mean estimate:
# 0 0 0 0 0 0 0 1 1 1 1



# Spectral Bias:
# 0.25 0.20 0.28 0.47 0.60 0.67 0.67 0.59 0.54 0.57




# Inputs :
# Run Time : 8.02804 mins
# n = 500
# iter = 1000
# R = 1
# B = 10
# Sim = 500
# Simulatenous Coverage (M Distance) :
# 0.832
# Simultaneous Coverage (Fisher's) ;
# 0.84
# Pointwise Coverage :
# 0.722 0.774 0.912 0.914 0.878 0.908 0.918 0.918 0.896 0.914
# Bias Beta :
# 0.450 0.282 0.300 0.380 0.406 0.452 0.470 0.488 0.492 0.466 0.480
# Spectral Bias :
# 0.262 0.254 0.384 0.534 0.590 0.612 0.650 0.594 0.588 0.570


# Inputs :
# Run Time : 9.848147 mins
# n = 500
# iter = 1000
# R = 5
# B = 10
# Sim = 500
# Simultaneous Coverage (M Distance):
# 0.674
# Simultaneous Coverage (Fisher's) ;
# 0.69
# Pointwise Coverage :
# 0.810 0.854 0.880 0.864 0.882 0.880 0.880 0.876 0.872 0.888
# Bias Beta :
# 0.460 0.406 0.380 0.404 0.442 0.476 0.464 0.456 0.514 0.520 0.504
# Spectral Bias :
# 0.366 0.392 0.524 0.540 0.544 0.546 0.526 0.564 0.510 0.576


# Inputs :
# Run Time : 5.176307 mins
# n = 500
# iter = 1000
# R = 5
# B = 5
# Sim = 500
# Simultaneous Coverage:
# 0.61
# Pointwise Convergage :
# 0.722 0.774 0.912 0.914 0.878 0.908 0.918 0.918 0.896 0.914


# Inputs :
# Run Time : 20.54113 mins (even with computer pausing for a minute or so)
# n = 1000
# iter = 1000
# R = 10
# B = 10
# Sim = 500
# Simultaneous Coverage:
# 0.464
# Pointwise Convergage :
# 0.812 0.842 0.840 0.852 0.838 0.844 0.852 0.848 0.866 0.818

# Try different using Fisher's F-ratio: Calculates where 95% of my points under
# this distance (scaled M-distance) will be closer than the value of qf below.
Frat = nrow(Beta) * (Beta_distance) * (nrow(Beta) - B - 1) / ((B+1) *(nrow(Beta) - 1) * (nrow(Beta) + 1))
critf = qf(0.95, B+1, nrow(Beta) - (B+1))
mean(Frat <= critf)
# 0.486

#######################################
# Diagnostic Plots
#######################################
mean(Simult_truespec)
# 53%

# Seed 100 :
# 0.723
# Functional Boxplots
# Create a functional boxplot for the estimated log of the spectral density
# Will show the most recent Simulation result
par(mfrow = c(1,1))
pdf("~/Documents/Tessellation/R/Simulation_Results/Credible_Simultaneous_Interval.pdf")
Simult_specdens = fbplot(fit = log(test_specdens), x = test_omega, xlim = range(test_omega), ylim = range(log(test_specdens)), prob = c(0.95, 0.5, .05),
       xlab = "omega",
       factor = 1.5,
       barcol = "dodgerblue" ,
       col = c("purple", "orange", "violetred3"),
       fullout = TRUE,
       main = "Credible Interval (Simultaneous)")


# select curves where 95% falls above this depth
deepcurve = log(test_specdens)[, Simult_specdens$depth >= quantile(Simult_specdens$depth, .05)]
uppercurve = apply(deepcurve, 1, FUN = max)
lowercurve = apply(deepcurve, 1, FUN = min)

# Functional boxplot 95% bounds
lines(test_omega, uppercurve, lwd = 2)
lines(test_omega, lowercurve, lwd = 2)
lines(test_omega, log(true_spec), lwd = 2, col = "green")
legend("topright",
       legend = c("middle 95", "middle 50", "middle 5", "true spectral density"),
       col = c("purple", "orange", "violetred3", "green"),
       pch = c(20,20,20,20))
dev.off()
# Calculate if truespec dens completely fall within uppercurve and lowercurve
# Will return TRUE if true curve is everywhere in bounds
all((log(true_spec) <= uppercurve) & (log(true_spec) >= lowercurve))
# seed31/ R1/ n = 500:
# TRUE
# seed100/ R1 / n = 500:
# TRUE

# seed 31 / R5 / n = 500 :
# FALSE
# seed 100 / R5 / n = 500 :
# TRUE

# seed 31 / R5 / n = 1000:
# FALSE
# seed 100 / R5 / n = 1000:
# TRUE

# Create Trace Plots using mcmc package
colnames(Theta) = c(paste("b", seq(0,B), sep = ""), "tausq")
# tau^2 trace plot
pdf(file = "~/Documents/Tessellation/R/Simulation_Results/tausq_traceplots.pdf")
mcmc_trace(Theta[,B+2, drop = FALSE])
dev.off()
pdf(file = "~/Documents/Tessellation/R/Simulation_Results/beta_traceplots.pdf")
mcmc_trace(Theta[,-c(B+2), drop = FALSE])
dev.off()


#############################################
# Plot the posterior predictive and true
#############################################
pdf(file = "~/Documents/Tessellation/R/Simulation_Results/posterior_predictive.pdf")
Simult_specdens = fbplot(fit = log(post_pred[,,Sim]), x = test_omega, xlim = range(test_omega), ylim = range(log(post_pred[,,Sim])), prob = c(0.95, 0.5, .05),
                         xlab = "omega",
                         factor = 1.5,
                         barcol = "dodgerblue" ,
                         col = c("purple", "orange", "violetred3"),
                         fullout = TRUE,
                         ylab = "periodogram",
                         main = "Posterior Predictive")
lines(x = omega, y = log(av_perio), col = "green", lwd = 2)
legend("topright",
       legend = c("middle 95", "middle 50", "middle 5", "true periodogram"),
       col = c("purple", "orange", "violetred3", "green"),
       pch = c(20,20,20,20))
dev.off()

#################################################
# mcmcse : not helpful unless altering iter count
#################################################

# tau^2
mcse(x = Theta[,B+2])
# est
# 1.0459
# se
# 0.0175

#beta^*

# not sure which function to use from mcse mat or multi
# roughly does standard error
mcse_beta = mcse.mat(x = Theta[,-c(B+2)])

round(cbind(mcse_beta, true_Beta), 5)
     # est           se
# 0.0171978037 0.0005397890
# 0.7198770277 0.0005203687
# 0.1712138572 0.0005795037
# 0.0795108797 0.0005537198
# 0.0325336251 0.0006688887
# -0.0007262838 0.0005110527
# -0.0104669329 0.0005844425
# 0.0045298895 0.0006378697
# -0.0199138846 0.0005719261
# 0.0048291834 0.0006735232
# 0.0019303400 0.0004736163

burnin

dim(tcrossprod(rep(1,iter - burnin),true_Beta))
# 90 11
dim(Theta[,-c(B+2)])
# 990 11
# RMSE
colMeans((Beta - tcrossprod(rep(1,iter - burnin),true_Beta))^2)


# plot(confRegion(mcse_beta, which = c(1,2), level = .95), type = 'l', asp = 1)
# lines(confRegion(mcse_beta, which = c(1,2), level = .90), col = "red")

# across all simulations how far away beta's are from "true values"
colMeans(rmse)
# n = 500
# iter = 1000
# R = 20
# B = 10
# phi = 0.5
# 0.01620574 0.01552675 0.01667892 0.01660985 0.01600331 0.01532764 0.01454676 0.01431904
# 0.01315319 0.01464295 0.01270879

