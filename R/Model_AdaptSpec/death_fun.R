# Death move function

# Function for Birth Move

source("~/Documents/Tessellation/R/Model_AdaptSpec/get_r_fun.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/dtau.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/rtau.R")
source("~/Documents/Tessellation/R/Model_AdaptSpec/betapar.R")
death_fun = function(xi_prop,xi_cur,tmin,Smax,Beta,tau, timeseries, sigmasalpha, D, B, nu_0){
  
  # q(S^p | S^c) * q(\xi^p | S^c,\xi^c)
  # just the probability birth is possible given by the function get_r_fun
  q_seg_xi_birth = get_r_fun(Smax = Smax, tmin = tmin, xi_cur = xi_cur, xi_prop = xi_prop)
  q_seg_xi_death = get_r_fun(Smax = Smax, tmin = tmin, xi_cur = xi_prop, xi_prop = xi_cur)
  
  # Determine q_tau
  # Propose a tau
  m_star = which(xi_cur != xi_prop[-length(xi_prop)])[1]-1
  tau_cur = tau[m_star]
  beta_cur = Beta[,m_star]
  # call function that proposes two new tau parameters
  tau_p1 = rtau(tau_cur)
  tau_p2 = tau_cur^2/tau_p1
  if(m_star == 1){
    tau_prop = c(tau_p1, tau_p2, tau[(m_star+1):length(tau)])
  } else if(m_star == length(tau)){
    tau_prop = c(tau[1:(m_star-1)], tau_p1, tau_p2)
  }else{
    tau_prop = c(tau[1:(m_star-1)], tau_p1, tau_p2, tau[(m_star+1):length(tau)])
  }
  # Calculate q_tau
  q_tau = dtau(tau_p1 = tau_p1, tau_cur = tau_cur)
  
  # Determine q_beta
  par1 = betapar(tsq = tau_p1, sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_prop[m_star]+1):(xi_prop[m_star+1])], B = B)
  par2 = betapar(tsq = tau_p2, sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_prop[m_star+1]+1):(xi_prop[m_star+2])], B = B)
  
  # get proposal betas
  beta_p1 = c(mvtnorm::rmvnorm(1, mean = par1$mean, sigma = par1$cov))
  beta_p2 = c(mvtnorm::rmvnorm(1, mean = par2$mean, sigma = par2$cov))
  
  # get q_betas
  q_beta1 = mvtnorm::dmvnorm(beta_p1, mean = par1$mean, sigma = par1$cov)
  q_beta2 = mvtnorm::dmvnorm(beta_p2, mean = par2$mean, sigma = par2$cov)
  
  # q_betac
  par3 = betapar(tsq = tau_cur, sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_prop[m_star]+1):(xi_prop[m_star+2])], B = B)
  q_betac = mvtnorm::dmvnorm(beta_cur, mean = par3$mean, sigma = par3$cov)
  
  
  # Calculate the Likelihood for birth
  # proposal:
  L_p1 = log_likelihood_adapt(b=beta_p1, Psi = par1$Psi, sumPsi = par1$sumPsi, perio = par1$perio) 
  L_p2 = log_likelihood_adapt(b=beta_p2, Psi = par2$Psi, sumPsi = par2$sumPsi, perio = par2$perio)
  
  # current:
  L_c = log_likelihood_adapt(b=beta_cur, Psi = par3$Psi, sumPsi = par3$sumPsi, perio = par3$perio)
  
  # Calculate posterior
  # proposal:
  
  post_prop = L_p1 + L_p2 +
    dxi(xi = xi_prop,tmin = tmin) +
    sum(dnorm(beta_p1,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_p1)), log = TRUE)) +
    sum(dnorm(beta_p2,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_p2)), log = TRUE)) +
    invgamma::dinvgamma(tau_p1, nu_0/2, nu_0, log = TRUE) +
    invgamma::dinvgamma(tau_p2, nu_0/2, nu_0, log = TRUE)
  
  post_cur = L_c +
    dxi(xi = xi_cur,tmin = tmin) +
    sum(dnorm(beta_cur,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_cur)), log = TRUE)) +
    invgamma::dinvgamma(tau_cur, nu_0/2, nu_0, log = TRUE)
  
  # Calculate acceptance ratio
  q_prop = q_seg_xi_birth + q_tau + q_beta1 + q_beta2
  q_cur = q_seg_xi_death + q_betac
  
  A = exp(post_prop + q_cur - (post_cur + q_prop))
  # Metropolis step:
  U = runif(1)
  if(U < A){
    # Accept
    if(m_star == 1){
      Beta_prop = cbind(beta_p1, beta_p2, Beta[(m_star+1):(ncol(Beta))])
    } else if(m_star == ncol(Beta)){
      Beta_prop = cbind(Beta[1:(m_star-1)], beta_p1, beta_p2)
    }else{
      Beta_prop = cbind(Beta[1:(m_star-1)], beta_p1, beta_p2, Beta[(m_star+1):(ncol(Beta))])
    }
    return(list("Beta" = Beta_prop, "xi" = xi_prop, "tau" = tau_prop))
  }else{
    # Reject
    return(list("Beta" = Beta, "xi" = xi_cur, "tau" = tau))
  }
  
}










