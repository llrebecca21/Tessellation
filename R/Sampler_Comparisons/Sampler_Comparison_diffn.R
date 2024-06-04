# Sampler Comparison different lengths of n for each replicate time series
library(RColorBrewer)
# Run a comparison between the different sampling algorithms and plot their mean spectral density on the same plot.
set.seed(105)
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
#phi = c(1.4256,
#        -0.7344,
#        0.1296)
R = 6
B = 10
avg_n = 1000


# Initialize a list to store each time series of length "n"
ts_list = vector(mode = "list", length = R)
for(r in 1:R){
  #n = rpois(n = 1, lambda = avg_n)
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
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(100,1 / (4 * pi * (1:B)^2))))

# plots of the covariance matrices
if(FALSE){
  # custom color function
  col1 = "blue"
  col2 = "yellow"
  n = 100
  col_fun = function(col1,col2,n){
    #create seq
    A = seq(0,1,len = n)
    A = tcrossprod(A,col2rgb(col1)[,1])+tcrossprod((1-A),col2rgb(col2)[,1])
    A = A/255
    apply(A,1,FUN = function(a){rgb(a[1],a[2],a[3])})
  }
  lighten = function(x, y = 0){
    
    x = col2rgb(x) / 255
    lo = log(x / (1 - x)) + y
    x = 1 / (1 + exp(-lo))
    rgb(x[1], x[2], x[3])
    
  }
  
  
  
  # Plots of the covariance matrix
  # eta_br
  gap = 4
  M = matrix(-Inf, 11, 11 + (R-1)*gap)
  for(i in 1:R){
    
    M[cbind(1:11, 1:11 + gap * (i - 1))] = c(1,Result_eta_br_n$Theta[1000,12]*Result_eta_br_n$D)/Result_eta_br_n$eta_br_array[1000,,i]
    
  }
  pdf(file = "eta_br_n_cov.pdf",
      width = 10,
      height = 10)
  par(mar=c(0,0,0,0))
  M = t(M[ , ])[ , 11:1]
  image(sign(M) * sqrt(abs(M)), col = sapply(brewer.pal(n = 9, name = "OrRd"), function(x) lighten(x, -1.5)),
        xaxs = NULL, yaxs = NULL, axes = FALSE)
  dev.off()
  
  #hcl.pals()
  
  # eta_r
  gap = 4
  M = matrix(-Inf, 11, 11 + (R - 1) * gap)
  for(i in 1:R){
    
    M[cbind(1:11, 1:11 + gap * (i - 1))] = c(1,Result_eta_r_n$Theta[1000,12]*Result_eta_r_n$D)/Result_eta_r_n$eta_array[1000,i]
    
  }
  pdf(file = "eta_r_n_cov.pdf",
      width = 10,
      height = 10)
  par(mar=c(0,0,0,0))
  M = t(M[ , ])[ , 11:1]
  image(sign(M) * sqrt(abs(M)), col = sapply(brewer.pal(n = 9, name = "OrRd"), function(x) lighten(x, -1.5)),
        xaxs = NULL, yaxs = NULL,axes = FALSE)
  dev.off()
  # hcl.colors
  # hcl.colors(100, palette = hcl.pals()[30],rev = TRUE)
  # col_fun("salmon","darkseagreen1",100)
  #dev.off()
  #21
  #24
  #30
  #44
  
  
  # Wishart
  M = solve(Result_Wishart_n$Lambda_array[1000,,])[,seq(B+1,1)]
  pdf(file = "Wishart_cov.pdf",
      width = 10,
      height = 10)
  par(mar=c(0,0,0,0))
  image(sign(M) * sqrt(abs(M)), col = sapply(brewer.pal(n = 9, name = "OrRd"), function(x) lighten(x, -1.5)),
        xaxs = NULL, yaxs = NULL, axes = FALSE)
  dev.off()
  
}

# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor(n / 2)
omega = (2 * pi * (1:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1


# true spectral density vs. models
true_spec = log(arma_spec(omega = omega, phi = phi))
# store values
eta_br_MSE = rep(NA,R)
eta_r_MSE = rep(NA,R)
wish_MSE = rep(NA,R)
ind_MSE = rep(NA,R)
com_MSE = rep(NA,R)
set.seed(101)
for(i in 1:R){
  # eta_br
  est_spec_etabr1 = (Psi %*% t(Result_eta_br_n$bb_beta_array[,,i]))
  eta_br_MSE[i] = mean(colMeans((est_spec_etabr1 - true_spec)^2))
  # eta_r
  est_spec_etar = (Psi %*% t(Result_eta_r_n$bb_beta_array[,,i]))
  eta_r_MSE[i] = mean(colMeans(est_spec_etar - true_spec)^2)
  # Wishart
  est_spec_wish = (Psi %*% t(Result_Wishart_n$bb_beta_array[,,i]))
  wish_MSE[i] = mean(colMeans(est_spec_wish - true_spec)^2)
  # Independent
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[i]]), B = B)
  est_spec_ind = (Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  ind_MSE[i] = mean(colMeans(est_spec_ind - true_spec)^2)
  # Common
  est_spec_common = (Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  com_MSE[i] =  mean(colMeans(est_spec_common - true_spec)^2)
}
# average across the replicates
mean(eta_br_MSE) # 0.0274
mean(eta_r_MSE) # 0.0094
mean(wish_MSE) # 0.0103
mean(ind_MSE) # 0.0079
mean(com_MSE) # 0.0003


# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
pdf(file = "Spec_Dens_AR2.pdf",
    width = 10,
    height = 10)
par(mfrow = c(3,1), mar = c(4,4,3,0.1)+0.1)
set.seed(105)
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
         lwd = c(2,1,1,1, 1), legend = c("True", "Common", "Independent","Unstructured", "Replicate","Replicate & Basis"))
}
dev.off()



# Plot the Posterior Mean Estimates with True Spectral Density
# cairo_pdf(file = "Post_Mean_AR2.pdf",
#     width = 10,
#     height = 7)
# set.seed(105)
# par(fig = c(0,0.4,.514,1), mai = c(0,1.1,0.8,0), las = 1)
# if(TRUE){
#   r = 1
#   specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
#   specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
#   Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
#   specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
#   specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
#   specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
#   plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-1.5,2), ylab = "log(Spectral Density)", xlab = expression(omega),
#        main = "Posterior Mean Estimates \nwith True Spectral Density", )
#   lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
#   lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
#   lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
#   lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
#   lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
#   lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2, lty = 2)
#   #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
#   legend("topright", col = c("black",
#                              rgb(.76, .65, .81),
#                              rgb(1, .412, .71),
#                              rgb(.48, .19, .58),
#                              rgb(.65, .85, .62),
#                              rgb(0, .53, .21)),
#          lwd = c(2,2,2,2,2,2),
#          lty = c(2,1,1,1,1,1),
#          legend = c("True", "Common", "Independent","Unstructured","Replicate","Replicate & Basis"),
#          bty = "n")
# }

cairo_pdf(file = "Post_Mean_AR2.pdf",
    width = 10,
    height = 4)
set.seed(105)
par(mfrow = c(1,3))
for(r in 1:3){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,2), ylab = "log(Spectral Density)", xlab = expression(omega),
       main = "Posterior Mean Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = phi)), col = "black", lwd = 2, lty = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  legend("topright", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(1, .412, .71),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(2,2,2,2,2,2),
         lty = c(2,1,1,1,1,1),
         legend = c("True", "Common", "Independent","Unstructured","Replicate","Replicate & Basis"),
         bty = "n")
}
dev.off()


###########################
# Run with generate_Krafty
###########################
R = 10
B = 10
set.seed(105)

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
J = floor(n / 2)
omega = (2 * pi * (1:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

#Plot the Spectral Density Estimates with True Spectral Density
pdf(file = "Spec_Dens_MA.pdf",
    width = 10,
    height = 10)
set.seed(105)
par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
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
dev.off()


# Plot the Posterior Mean Estimates with True Spectral Density
pdf(file = "Post_Mean_MA.pdf",
    width = 10,
    height = 10)
set.seed(105)
par(mfrow = c(2,2))
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
dev.off()

# true spectral density vs. models
# store values
eta_br_MSE = rep(NA,R)
eta_r_MSE = rep(NA,R)
wish_MSE = rep(NA,R)
ind_MSE = rep(NA,R)
com_MSE = rep(NA,R)
set.seed(105)
for(i in 1:R){
  # true spec
  #true_spec = log(arma_spec(omega = omega,
  #                          phi = gendata1[[i]]$phi) + arma_spec(omega = omega, phi = gendata2[[i]]$phi))
  true_spec = log(arma_spec(omega = omega, theta = theta_true[i]))
  # eta_br
  est_spec_etabr1 = (Psi %*% t(Result_eta_br_n$bb_beta_array[,,i]))
  eta_br_MSE[i] = mean(colMeans((est_spec_etabr1 - true_spec)^2))
  # eta_r
  est_spec_etar = (Psi %*% t(Result_eta_r_n$bb_beta_array[,,i]))
  eta_r_MSE[i] = mean(colMeans(est_spec_etar - true_spec)^2)
  # Wishart
  est_spec_wish = (Psi %*% t(Result_Wishart_n$bb_beta_array[,,i]))
  wish_MSE[i] = mean(colMeans(est_spec_wish - true_spec)^2)
  # Independent
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[i]]), B = B)
  est_spec_ind = (Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  ind_MSE[i] = mean(colMeans(est_spec_ind - true_spec)^2)
  # Common
  est_spec_common = (Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  com_MSE[i] =  mean(colMeans(est_spec_common - true_spec)^2)
}
# average across the replicates
mean(eta_br_MSE) # 0.8016
mean(eta_r_MSE) # 0.0314
mean(wish_MSE) # 0.0067
mean(ind_MSE) # 0.0080
mean(com_MSE) # 1.8291



###########################
# Generate an AR(2) Mix
###########################
set.seed(105)
R = 6
B = 10
#avg_n = 750

peaks1 <- runif(R, min = 0.19, max = 0.21)
bandwidths1 <- rep(0.1, R)
peaks2 <- runif(R, min = (pi / 4) - 0.01, max = (pi/4) + 0.01)
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
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(100,1 / (4 * pi * (1:B)^2))))


# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor(n / 2)
omega = (2 * pi * (1:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density

set.seed(105)
par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
for(r in 1:3){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,10), ylab = "Spectral Density", xlab = expression(omega),
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
  # legend("topright", col = c("black",
  #                            #Single light purple
  #                            rgb(.76, .65, .81),
  #                            #Independent hot pink
  #                            rgb(1, .412, .71, 0.4),
  #                            #Wishart dark purple
  #                            rgb(.48, .19, .58),
  #                            #eta_r light green
  #                            rgb(.65, .85, .62),
  #                            #eta_br dark green
  #                            rgb(0, .53, .21)),
  #        lwd = c(2,1,1,1, 1), legend = c("True", "Common", "Independent","Unstructured", "Replication","Replication and Basis"))
}



# Plot the Posterior Mean Estimates with True Spectral Density
cairo_pdf(file = "Post_mean_ar2_boring.pdf",
          width = 10,
          height = 4)
par(mfrow = c(1,3))
set.seed(105)
for(r in 1:3){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,8), ylab = "log(Spectral Density)", xlab = expression(omega),
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "gray", lwd = 0.5)
  # legend("topright", col = c("black",
  #                            rgb(.76, .65, .81),
  #                            rgb(1, .412, .71),
  #                            rgb(.48, .19, .58),
  #                            rgb(.65, .85, .62),
  #                            rgb(0, .53, .21)),
  #        lwd = c(3,3,3,3,3),
  #        legend = c("True", "Single", "Independent","Wishart","eta_r","eta_br"))
}
dev.off()

# true spectral density vs. models
# store values
eta_br_MSE = rep(NA,R)
eta_r_MSE = rep(NA,R)
wish_MSE = rep(NA,R)
ind_MSE = rep(NA,R)
com_MSE = rep(NA,R)
set.seed(105)
for(i in 1:R){
  # true spec
  true_spec = log(arma_spec(omega = omega,
                            phi = gendata1[[i]]$phi) + arma_spec(omega = omega, phi = gendata2[[i]]$phi))
  # eta_br
  est_spec_etabr1 = (Psi %*% t(Result_eta_br_n$bb_beta_array[,,i]))
  eta_br_MSE[i] = mean(colMeans((est_spec_etabr1 - true_spec)^2))
  # eta_r
  est_spec_etar = (Psi %*% t(Result_eta_r_n$bb_beta_array[,,i]))
  eta_r_MSE[i] = mean(colMeans(est_spec_etar - true_spec)^2)
  # Wishart
  est_spec_wish = (Psi %*% t(Result_Wishart_n$bb_beta_array[,,i]))
  wish_MSE[i] = mean(colMeans(est_spec_wish - true_spec)^2)
  # Independent
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[i]]), B = B)
  est_spec_ind = (Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  ind_MSE[i] = mean(colMeans(est_spec_ind - true_spec)^2)
  # Common
  est_spec_common = (Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  com_MSE[i] =  mean(colMeans(est_spec_common - true_spec)^2)
}
# average across the replicates
mean(eta_br_MSE) # 0.1395
mean(eta_r_MSE) # 0.0311
mean(wish_MSE) # 0.0256
mean(ind_MSE) # 0.0336
mean(com_MSE) # 0.0233

# Mix_AR2_2 ---------------------------------------------------------------

###########################
# Generate an AR(2) Mix 
# Example 2
###########################
set.seed(105)
R = 20
B = 25
#avg_n = 750

peaks1 <- runif(R, min = 2.3, max = 2.5)
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


par(mfrow = c(2,3))
plot(x = omega, y = log(arma_spec(omega = omega,
                                  phi = gendata1[[1]]$phi) + arma_spec(omega = omega, phi = gendata2[[1]]$phi)), col = "black", lwd = 2)

plot(x = omega, y = log(arma_spec(omega = omega,
                                  phi = gendata1[[2]]$phi) + arma_spec(omega = omega, phi = gendata2[[2]]$phi)), col = "black", lwd = 2)

plot(x = omega, y = log(arma_spec(omega = omega,
                                  phi = gendata1[[3]]$phi) + arma_spec(omega = omega, phi = gendata2[[3]]$phi)), col = "black", lwd = 2)
plot(x = omega, y = log(arma_spec(omega = omega,
                                  phi = gendata1[[4]]$phi) + arma_spec(omega = omega, phi = gendata2[[4]]$phi)), col = "black", lwd = 2)


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
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(100,1 / (4 * pi * (1:B)^2))))


# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor(n / 2)
omega = (2 * pi * (1:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# Plot true Spec Density
plot(x = omega, y = log(arma_spec(omega = omega,
                                   phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
set.seed(105)
par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
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
  lines(x = omega, y = log(arma_spec(omega = omega,
                                     phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
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
set.seed(105)
par(mfrow = c(2,2))
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

# true spectral density vs. models
# store values
eta_br_MSE = rep(NA,R)
eta_r_MSE = rep(NA,R)
wish_MSE = rep(NA,R)
ind_MSE = rep(NA,R)
com_MSE = rep(NA,R)
set.seed(105)
for(i in 1:R){
  # true spec
  true_spec = log(arma_spec(omega = omega,
                            phi = gendata1[[i]]$phi) + arma_spec(omega = omega, phi = gendata2[[i]]$phi))
  # eta_br
  est_spec_etabr1 = (Psi %*% t(Result_eta_br_n$bb_beta_array[,,i]))
  eta_br_MSE[i] = mean(colMeans((est_spec_etabr1 - true_spec)^2))
  # eta_r
  est_spec_etar = (Psi %*% t(Result_eta_r_n$bb_beta_array[,,i]))
  eta_r_MSE[i] = mean(colMeans(est_spec_etar - true_spec)^2)
  # Wishart
  est_spec_wish = (Psi %*% t(Result_Wishart_n$bb_beta_array[,,i]))
  wish_MSE[i] = mean(colMeans(est_spec_wish - true_spec)^2)
  # Independent
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[i]]), B = B)
  est_spec_ind = (Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  ind_MSE[i] = mean(colMeans(est_spec_ind - true_spec)^2)
  # Common
  est_spec_common = (Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  com_MSE[i] =  mean(colMeans(est_spec_common - true_spec)^2)
}
# average across the replicates
mean(eta_br_MSE) # 0.8016
mean(eta_r_MSE) # 0.0314
mean(wish_MSE) # 0.0067
mean(ind_MSE) # 0.0080
mean(com_MSE) # 1.8291






# Mix_AR3 -----------------------------------------------------------------


###########################
# Generate an AR(2) Mix 
# Example 3
###########################
set.seed(105)
R = 4
B = 25
#avg_n = 750

peaks1 <- runif(R, min = 0.5, max = 1)
bandwidths1 <- rep(0.05, R)
# peaks2 <- runif(R, min = (pi / 4) - 0.2, max = (pi/4) + 0.2)
# bandwidths2 <- rep(0.1, R)
peaks2 <- runif(R, min = 2, max = 2.5)
bandwidths2 <- rep(0.03, R)


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

# Plot true Spec Density
plot(x = omega, y = log(arma_spec(omega = omega,
                                  phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))


#Plot the Spectral Density Estimates with True Spectral Density
par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
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
  lines(x = omega, y = log(arma_spec(omega = omega,
                                     phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topleft", col = c("black",
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
par(mfrow = c(2,2))
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
  legend("topleft", col = c("black",
                             rgb(.76, .65, .81),
                             rgb(1, .412, .71),
                             rgb(.48, .19, .58),
                             rgb(.65, .85, .62),
                             rgb(0, .53, .21)),
         lwd = c(3,3,3,3,3),
         legend = c("True", "Single", "Independent","Wishart","eta_r","eta_br"))
}



#Plot the Spectral Density Estimates with True Spectral Density : nonlog
par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
for(r in 1:R){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(0,1000), ylab = "Spectral Density", xlab = "omega",
       main = "Spectral Density Estimates \nwith True Spectral Density")
  # Plot Model Single n
  for(h in sample(ncol(specdens_Single_n), 100, replace = FALSE)){ #light purple
    lines(x = omega, y = (specdens_Single_n[,h]), col = rgb(.76, .65, .81, 0.4))
  }
  # Plot Independent Model n
  for(h in sample(ncol(specdens_Independent_n), 100, replace = FALSE)){ #hotpink
    lines(x = omega, y = (specdens_Independent_n[,h]), col = rgb(1, .412, .71, 0.4))
  }
  # Plot Model Wishart n
  for(h in sample(ncol(specdens_Wishart_n), 100, replace = FALSE)){ #dark purple
    lines(x = omega, y = (specdens_Wishart_n[,h]), col = rgb(.48, .19, .58, 0.4))
  }
  # Plot Model eta_r
  for(h in sample(ncol(specdens_eta_r_n), 100, replace = FALSE)){ #light green
    lines(x = omega, y = (specdens_eta_r_n[,h]), col = rgb(.65, .85, .62, 0.4))
  }
  # Plot Model eta_br_n
  for(h in sample(ncol(specdens_eta_br_n), 100, replace = FALSE)){ #dark green
    lines(x = omega, y = (specdens_eta_br_n[,h]), col = rgb(0, .53, .21, 0.4))
  }
  lines(x = omega, y = (arma_spec(omega = omega,
                                     phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
  #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
  legend("topleft", col = c("black",
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




# AR2_4 -------------------------------------------------------------------

###########################
# Generate an AR(2) Mix 
# Example 4
###########################
set.seed(105)
R = 6
B = 10
#avg_n = 750

peaks1 <- runif(R, min = 0.2, max = 2.5)
bandwidths1 <- rep(0.05, R)
# peaks2 <- runif(R, min = (pi / 4) - 0.2, max = (pi/4) + 0.2)
# bandwidths2 <- rep(0.1, R)
peaks2 <- runif(R, min = 2.0, max = 3.0)
bandwidths2 <- rep(1, R)


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

# gendata1: AR(2) coefficients for 'first' "peak" in each time series
# gendata2: AR(2) coefficients for 'second' "peak" in each time series
# gendata1[[1]] # phi = c(1.7330533,-0.9048374)
# gendata2[[1]] # phi = c(-0.6120395, -0.1353353)
# 
# gendata1[[2]] # phi = c(-1.5126316, -0.9048374)
# gendata2[[2]] # phi = c(-0.6921953, -0.1353353)
# 
# gendata1[[3]] # phi = c(1.1262363,-0.9048374)
# gendata2[[3]] # phi = c(-0.7147854, -0.1353353)
# 
# gendata1[[4]] # phi = c(0.8378112, -0.9048374)
# gendata2[[4]] # phi = c(-0.6998616, -0.1353353)

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
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(100,1 / (4 * pi * (1:B)^2))))


# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor(n/2)
omega = (2 * pi * (1:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

# Plot true Spec Density
par(mfrow = c(2,3))
plot(x = omega, y = log(arma_spec(omega = omega,
                                    phi = gendata1[[1]]$phi) + arma_spec(omega = omega, phi = gendata2[[1]]$phi)), col = "black", lwd = 2)

plot(x = omega, y = log(arma_spec(omega = omega,
                                    phi = gendata1[[2]]$phi) + arma_spec(omega = omega, phi = gendata2[[2]]$phi)), col = "black", lwd = 2)
  
plot(x = omega, y = log(arma_spec(omega = omega,
                                    phi = gendata1[[3]]$phi) + arma_spec(omega = omega, phi = gendata2[[3]]$phi)), col = "black", lwd = 2)
plot(x = omega, y = log(arma_spec(omega = omega,
                                    phi = gendata1[[4]]$phi) + arma_spec(omega = omega, phi = gendata2[[4]]$phi)), col = "black", lwd = 2)
  

# solve(crossprod(Psi), crossprod(Psi, log(Result_Single$perio)))
  
# true spectral density vs. models
# store values
eta_br_MSE = rep(NA,R)
eta_r_MSE = rep(NA,R)
wish_MSE = rep(NA,R)
ind_MSE = rep(NA,R)
com_MSE = rep(NA,R)
set.seed(105)
for(i in 1:R){
  # true spec
  true_spec = log(arma_spec(omega = omega,
                            phi = gendata1[[i]]$phi) + arma_spec(omega = omega, phi = gendata2[[i]]$phi))
  # eta_br
  est_spec_etabr1 = (Psi %*% t(Result_eta_br_n$bb_beta_array[,,i]))
  eta_br_MSE[i] = mean(colMeans((est_spec_etabr1 - true_spec)^2))
  # eta_r
  est_spec_etar = (Psi %*% t(Result_eta_r_n$bb_beta_array[,,i]))
  eta_r_MSE[i] = mean(colMeans(est_spec_etar - true_spec)^2)
  # Wishart
  est_spec_wish = (Psi %*% t(Result_Wishart_n$bb_beta_array[,,i]))
  wish_MSE[i] = mean(colMeans(est_spec_wish - true_spec)^2)
  # Independent
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[i]]), B = B)
  est_spec_ind = (Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  ind_MSE[i] = mean(colMeans(est_spec_ind - true_spec)^2)
  # Common
  est_spec_common = (Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  com_MSE[i] =  mean(colMeans(est_spec_common - true_spec)^2)
}
# average across the replicates
mean(eta_br_MSE) # 0.8016
mean(eta_r_MSE) # 0.0314
mean(wish_MSE) # 0.0067
mean(ind_MSE) # 0.0080
mean(com_MSE) # 1.8291


#Plot the Spectral Density Estimates with True Spectral Density
pdf(file = "Spec_Dens_weird.pdf",
    width = 10,
    height = 10)
par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
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
  lines(x = omega, y = log(arma_spec(omega = omega,
                                     phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
  J_perio = length(c(Result_Wishart_n$perio[[r]]))
  omega_perio = (2 * pi * (1:J_perio)) / list_n[r]
  points(x = omega_perio, y = log(Result_Wishart_n$perio[[r]]), col = "green", lwd = 0.5)
  legend("topleft", col = c("black",
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
         lwd = c(2,1,1,1, 1), legend = c("True", "Common", "Independent","Unstructured", "Replicate","Replicate & Basis"))
}
dev.off()

# Plot the Posterior Mean Estimates with True Spectral Density
cairo_pdf(file = "Post_mean_weird.pdf",
    width = 10,
    height = 4)
set.seed(105)
par(mfrow = c(1,3))
for(r in 1:3){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-2,8), ylab = "log(Spectral Density)", xlab = expression(omega),
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = log(arma_spec(omega = omega, phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2, lty = 2)
  # plot the periodogram
  # J_perio = length(c(Result_Wishart_n$perio[[r]]))-1
  # omega_perio = (2 * pi * (0:J_perio)) / list_n[r]
  # points(x = omega_perio, y = log(Result_Wishart_n$perio[[r]]), col = "green", lwd = 0.5)
  # legend("topleft", col = c("black",
  #                           rgb(.76, .65, .81),
  #                           rgb(1, .412, .71),
  #                           rgb(.48, .19, .58),
  #                           rgb(.65, .85, .62),
  #                           rgb(0, .53, .21)),
  #        lwd = c(2,2,2,2,2,2),
  #        lty = c(2,1,1,1,1,1),
  #        legend = c("True", "Common", "Independent","Unstructured", "Replicate","Replicate & Basis"))
}
dev.off()

if(TRUE){
  #Plot the Spectral Density Estimates with True Spectral Density : nonlog
  pdf(file = "Spec_Dens_weird_woutlog.pdf",
      width = 10,
      height = 10)
  par(mfrow = c(2,2), mar = c(4,4,3,0.1)+0.1)
  for(r in 1:R){
    specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
    specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
    Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
    specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
    specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
    specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
    plot(x =c(), y=c(), xlim = c(0,3), ylim = c(0,1000), ylab = "Spectral Density", xlab = "omega",
         main = "Spectral Density Estimates \nwith True Spectral Density")
    # Plot Model Single n
    for(h in sample(ncol(specdens_Single_n), 100, replace = FALSE)){ #light purple
      lines(x = omega, y = (specdens_Single_n[,h]), col = rgb(.76, .65, .81, 0.4))
    }
    # Plot Independent Model n
    for(h in sample(ncol(specdens_Independent_n), 100, replace = FALSE)){ #hotpink
      lines(x = omega, y = (specdens_Independent_n[,h]), col = rgb(1, .412, .71, 0.4))
    }
    # Plot Model Wishart n
    for(h in sample(ncol(specdens_Wishart_n), 100, replace = FALSE)){ #dark purple
      lines(x = omega, y = (specdens_Wishart_n[,h]), col = rgb(.48, .19, .58, 0.4))
    }
    # Plot Model eta_r
    for(h in sample(ncol(specdens_eta_r_n), 100, replace = FALSE)){ #light green
      lines(x = omega, y = (specdens_eta_r_n[,h]), col = rgb(.65, .85, .62, 0.4))
    }
    # Plot Model eta_br_n
    for(h in sample(ncol(specdens_eta_br_n), 100, replace = FALSE)){ #dark green
      lines(x = omega, y = (specdens_eta_br_n[,h]), col = rgb(0, .53, .21, 0.4))
    }
    lines(x = omega, y = (arma_spec(omega = omega,
                                    phi = gendata1[[r]]$phi) + arma_spec(omega = omega, phi = gendata2[[r]]$phi)), col = "black", lwd = 2)
    #points(x = omega, y = log(Result_Wishart$perio[,r]), col = "green", lwd = 0.5)
    legend("topleft", col = c("black",
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
  dev.off()
}

# check the eta_r and eta_br in this simulation

colMeans(Result_eta_r_n$eta_array)

dim(Result_eta_br_n$eta_br_array)

eta_br_w = apply(Result_eta_br_n$eta_br_array, c(2,3), mean)

plot(log(Result_eta_br_n$eta_br_array[,2,1]), type = "l")


dim(Result_eta_br_n$Theta)

plot(Result_eta_br_n$Theta[,27], type = "l")

# identifiability issue
Result_eta_br_n$eta_br_array[,-1,1]/Result_eta_br_n$Theta[,27]
plot((Result_eta_br_n$eta_br_array[,-1,1]/Result_eta_br_n$Theta[,27])[,1], type ="l")

y_max = apply(Result_eta_br_n$eta_br_array[,-1,]/Result_eta_br_n$Theta[,27], 2, max)

# First time series
par(mfrow = c(5,5))
for(i in 1:25){
  plot((Result_eta_br_n$eta_br_array[,-1,1]/Result_eta_br_n$Theta[,27])[,i], type ="l", ylab = "eta_br/tausw",
  ylim = c(0,y_max[i]))
}

# Second Time Series
par(mfrow = c(5,5))
for(i in 1:25){
  plot((Result_eta_br_n$eta_br_array[,-1,2]/Result_eta_br_n$Theta[,27])[,i], type ="l", ylab = "eta_br/tausw",
  ylim = c(0,y_max[i]))
}

# Third Time Series
par(mfrow = c(5,5))
for(i in 1:25){
  plot((Result_eta_br_n$eta_br_array[,-1,3]/Result_eta_br_n$Theta[,27])[,i], type ="l", ylab = "eta_br/tausw",
       ylim = c(0,y_max[i]))
}

plot(Result_eta_br_n$bb_beta_array[,2,1], type = "l")

# Beta trace plots
# First time series
par(mfrow = c(5,5))
for(i in 1:25){
  plot(Result_eta_br_n$bb_beta_array[,i,1], type ="l", ylab = "beta")
  lines(Result_eta_r_n$bb_beta_array[,i,1], type = "l", col = "red")
}

# Second time series
par(mfrow = c(5,5))
for(i in 1:25){
  plot(Result_eta_br_n$bb_beta_array[,i,2], type ="l", ylab = "beta")
  lines(Result_eta_r_n$bb_beta_array[,i,2], type = "l", col = "red")
}


# Third time series
par(mfrow = c(5,5))
for(i in 0:25){
  plot(Result_eta_br_n$bb_beta_array[,i,3], type ="l", ylab = "beta")
  lines(Result_eta_r_n$bb_beta_array[,i,3], type = "l", col = "red")
}


# beta_generation_example -------------------------------------------------
set.seed(105)
R = 6
mother_beta = rep(0,11)
mother_beta[4] = 2

mother_beta_sd = rep(0.5,11)
mother_beta_sd[4] = 0.08
mother_beta_sd[1] = 0.08
# create true betas
beta_true = rmvnorm(R, mean = mother_beta, sigma = diag(mother_beta_sd^2)+tcrossprod(mother_beta_sd))
# plot generates true spectral densities from some mother beta
for(r in 1:R){
  plot(omega, Psi %*% beta_true[r,])
}

# generate time series
ts_list = vector(mode = "list", length = R)

# generate values for length of each time series
for(r in 1:R){
  n = sample(500:1000, 1)
  ts_list[[r]] = t(rX_from_specdens(n=1,n_times = n, basis = beta_true[r,]))
}

# plot the time series generated:
for(r in 1:R){
  ts.plot(ts_list[[r]])
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
Result_Wishart_n = Sampler_Wishart_n(ts_list = ts_list, B = B, tausquared = 1, V = diag(10/c(100,1 / (4 * pi * (1:B)^2))))


# Plot the results below
# Define omega and Psi for plotting purposes
n = 1000
J = floor(n/2)
omega = (2 * pi * (1:J)) / n
Psi = outer(X = omega, Y = 0:B, FUN = function(x,y){sqrt(2)* cos(y * x)})
# redefine the first column to be 1's
Psi[,1] = 1

cairo_pdf(file = "Post_mean_beta_gen.pdf",
          width = 10,
          height = 4)
set.seed(105)
par(mfrow = c(1,3))
for(r in 1:3){
  specdens_eta_br_n = exp(Psi %*% t(Result_eta_br_n$bb_beta_array[,,r]))
  specdens_eta_r_n = exp(Psi %*% t(Result_eta_r_n$bb_beta_array[,,r]))
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[r]]), B = B)
  specdens_Independent_n = exp(Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  specdens_Single_n = exp(Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  specdens_Wishart_n = exp(Psi %*% t(Result_Wishart_n$bb_beta_array[,,r]))
  plot(x =c(), y=c(), xlim = c(0,3), ylim = c(-5,8), ylab = "log(Spectral Density)", xlab = expression(omega),
       main = "Spectral Density Estimates \nwith True Spectral Density")
  lines(x = omega, y = rowMeans(log(specdens_Single_n)), lwd = 3, col = rgb(.76, .65, .81))
  lines(x = omega, y = rowMeans(log(specdens_Independent_n)), lwd = 3, col = rgb(1, .412, .71))
  lines(x = omega, y = rowMeans(log(specdens_Wishart_n)), lwd = 3, col = rgb(.48, .19, .58))
  lines(x = omega, y = rowMeans(log(specdens_eta_r_n)), lwd = 3, col = rgb(.65, .85, .62))
  lines(x = omega, y = rowMeans(log(specdens_eta_br_n)), lwd = 3, col = rgb(0, .53, .21))
  lines(x = omega, y = Psi %*% beta_true[r,], col = "black", lwd = 2, lty = 2)
  # plot the periodogram
  # J_perio = length(c(Result_Wishart_n$perio[[r]]))-1
  # omega_perio = (2 * pi * (0:J_perio)) / list_n[r]
  # points(x = omega_perio, y = log(Result_Wishart_n$perio[[r]]), col = "green", lwd = 0.5)
  # legend("topleft", col = c("black",
  #                           rgb(.76, .65, .81),
  #                           rgb(1, .412, .71),
  #                           rgb(.48, .19, .58),
  #                           rgb(.65, .85, .62),
  #                           rgb(0, .53, .21)),
  #        lwd = c(2,2,2,2,2,2),
  #        lty = c(2,1,1,1,1,1),
  #        legend = c("True", "Common", "Independent","Unstructured", "Replicate","Replicate & Basis"))
}
dev.off()



# true spectral density vs. models
# store values
eta_br_MSE = rep(NA,R)
eta_r_MSE = rep(NA,R)
wish_MSE = rep(NA,R)
ind_MSE = rep(NA,R)
com_MSE = rep(NA,R)
set.seed(105)
for(i in 1:R){
  # true spec
  true_spec = c(Psi %*% beta_true[i,])
  # eta_br
  est_spec_etabr1 = (Psi %*% t(Result_eta_br_n$bb_beta_array[,,i]))
  eta_br_MSE[i] = mean(colMeans((est_spec_etabr1 - true_spec)^2))
  # eta_r
  est_spec_etar = (Psi %*% t(Result_eta_r_n$bb_beta_array[,,i]))
  eta_r_MSE[i] = mean(colMeans(est_spec_etar - true_spec)^2)
  # Wishart
  est_spec_wish = (Psi %*% t(Result_Wishart_n$bb_beta_array[,,i]))
  wish_MSE[i] = mean(colMeans(est_spec_wish - true_spec)^2)
  # Independent
  Result_Independent_n = Sampler_Single_n(ts_list = list(ts_list[[i]]), B = B)
  est_spec_ind = (Psi %*% t(Result_Independent_n$Theta[,-(B+2)]))
  ind_MSE[i] = mean(colMeans(est_spec_ind - true_spec)^2)
  # Common
  est_spec_common = (Psi %*% t(Result_Single_n$Theta[,-(B+2)]))
  com_MSE[i] =  mean(colMeans(est_spec_common - true_spec)^2)
}
# average across the replicates
mean(eta_br_MSE) # 1.3459
mean(eta_r_MSE) # 0.1017
mean(wish_MSE) # 0.0922
mean(ind_MSE) # 0.1078
mean(com_MSE) # 1.1905



