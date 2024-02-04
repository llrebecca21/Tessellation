# Sampler Comparison different lengths of n for each replicate non-stationary time series

# Run a comparison between the different sampling algorithms and plot their mean spectral density on the same plot.
set.seed(105)
par(mfrow = c(2,4))
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

# First run with AR(p) data generating function as the input for both Samplers
phi1 = c(0.5, -0.4)
phi2 = c(0.8, -0.2)
#phi = 0.5
R = 1
B = 10
avg_n = 750


# Initialize a list to store each time series of length "n"
ts_list = vector(mode = "list", length = R)
for(r in 1:R){
  # n = sample(500:1000, 1)
  n = avg_n
  ts_list[[r]] = generate_nonstat_abrupt(phi1 = phi1, phi2 = phi2 ,n = n, R = 1)$matrix_timeseries
}

# plot the resulting time series
ts.plot(ts_list[[1]])






