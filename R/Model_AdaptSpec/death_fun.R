# Death move function

# Function for Birth Move

source("~/Documents/Tessellation/R/Model_AdaptSpec/get_r_fun.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/dtau.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/rtau.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/betapar.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/rxi_death.R")
death_fun = function(xi_cur,tmin,Smax,Beta,tau, timeseries, sigmasalpha, D, B, nu_0){
  # |xi_prop| < |xi_cur|
  # get proposal xi
  xi_prop = rxi_death(xi_cur = xi_cur)
  # q(S^p | S^c) * q(\xi^p | S^c,\xi^c)
  # just the probability death is possible given by the function get_r_fun
  q_seg_xi_death = get_r_fun(Smax = Smax, tmin = tmin, xi_cur = xi_cur, xi_prop = xi_prop)
  q_seg_xi_birth = get_r_fun(Smax = Smax, tmin = tmin, xi_cur = xi_prop, xi_prop = xi_cur)
  
  # Determine q_tau
  # Propose a tau
  # xi_cur = seq(0,1000,by = 100)
  # xi_prop = c(xi_cur[1:5],xi_cur[7:11])
  m_star = which(xi_prop != xi_cur[-length(xi_cur)])[1]-1
  
  # extract tau_c1
  tau_c1 = tau[m_star]
  # extract tau_c2
  tau_c2 = tau[m_star + 1]
  tau_p = sqrt(tau_c1 * tau_c2)
  
  if(length(xi_cur) == 3){
    tau_prop = tau_p
  }else if(m_star == 1){
    tau_prop = c(tau_p, tau[(m_star+2):length(tau)])
  } else if(m_star == (length(tau)-1)){
    tau_prop = c(tau[1:(m_star-1)], tau_p)
  }else{
    tau_prop = c(tau[1:(m_star-1)], tau_p, tau[(m_star+2):length(tau)])
  }
  
  q_tau = dtau(tau_p1 = tau_c1, tau_cur = tau_p)
  
  # extract beta_c1
  beta_c1 = Beta[m_star, ]
  # extract beta_c2
  beta_c2 = Beta[m_star + 1, ]

  # Determine q_beta
  par_c1 = betapar(tsq = tau_c1, sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_cur[m_star]+1):(xi_cur[m_star + 1])], B = B)
  par_c2 = betapar(tsq = tau_c2, sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_cur[m_star+1]+1):(xi_cur[m_star+2])], B = B)
  
  # get q_betas
  q_beta1 = mvtnorm::dmvnorm(beta_c1, mean = par_c1$mean, sigma = par_c1$cov, log = TRUE)
  q_beta2 = mvtnorm::dmvnorm(beta_c2, mean = par_c2$mean, sigma = par_c2$cov, log = TRUE)
  
  # q_betap
  par_p = betapar(tsq = tau_p, sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_cur[m_star]+1):(xi_cur[m_star+2])], B = B)
  beta_p = c(mvtnorm::rmvnorm(1, mean = par_p$mean, sigma = par_p$cov))
  q_betap = mvtnorm::dmvnorm(beta_p, mean = par_p$mean, sigma = par_p$cov, log = TRUE)
  
  # Calculate the Likelihood for birth
  # current:
  L_c1 = log_likelihood_adapt(b=beta_c1, Psi = par_c1$Psi, sumPsi = par_c1$sumPsi, perio = par_c1$perio)
  L_c2 = log_likelihood_adapt(b=beta_c2, Psi = par_c2$Psi, sumPsi = par_c2$sumPsi, perio = par_c2$perio)
  
  # proposal:
  L_p = log_likelihood_adapt(b=beta_p, Psi = par_p$Psi, sumPsi = par_p$sumPsi, perio = par_p$perio)
  
  # Calculate posterior
  # current:
  
  post_cur = L_c1 + L_c2 +
    dxi(xi = xi_cur,tmin = tmin) +
    sum(dnorm(beta_c1,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_c1)), log = TRUE)) +
    sum(dnorm(beta_c2,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_c2)), log = TRUE)) +
    invgamma::dinvgamma(tau_c1, nu_0/2, nu_0, log = TRUE) +
    invgamma::dinvgamma(tau_c2, nu_0/2, nu_0, log = TRUE)
  
  post_prop = L_p +
    dxi(xi = xi_prop,tmin = tmin) +
    sum(dnorm(beta_p,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_p)), log = TRUE)) +
    invgamma::dinvgamma(tau_p, nu_0/2, nu_0, log = TRUE)
  
  # Calculate acceptance ratio
  q_cur = q_seg_xi_birth + q_tau + q_beta1 + q_beta2
  q_prop = q_seg_xi_death + q_betap
  
  A = 1/exp(post_prop + q_cur - (post_cur + q_prop))
  # Metropolis step:
  U = runif(1)
  if(U < A){
    # Accept
    if(length(xi_cur) == 3){
      Beta_prop = matrix(beta_p, ncol = length(beta_p), nrow = 1)
    }else if(m_star == 1){
      Beta_prop = rbind(beta_p, Beta[(m_star+2):nrow(Beta), ])
    } else if(m_star == (nrow(Beta)-1)){
      Beta_prop = rbind(Beta[1:(m_star-1), ], beta_p)
    }else{
      Beta_prop = rbind(Beta[1:(m_star-1), ], beta_p, Beta[(m_star+2):nrow(Beta), ])
    }
    return(list("Beta" = Beta_prop, "xi" = xi_prop, "tau" = tau_prop))
  }else{
    # Reject
    return(list("Beta" = Beta, "xi" = xi_cur, "tau" = tau))
  }
  
}










