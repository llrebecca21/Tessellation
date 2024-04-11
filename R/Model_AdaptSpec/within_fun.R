# Function to perform a within step
source("~/Documents/Tessellation/R/Model_AdaptSpec/dxi_within.R")
within_fun = function(xi_cur,tmin,Smax,Beta,tau, timeseries, sigmasalpha, D, B, nu_0){
  # get the proposed xi_prop to move
  output = rxi_within(tmin = tmin, xi_cur = xi_cur)
  xi_prop = output$xi_prop
  m_star = output$m_star
  
  # get q_xi's
  q_xip = dxi_within(tmin = tmin, xi_cur = xi_cur, xi_prop = xi_prop)
  q_xic = dxi_within(tmin = tmin, xi_cur = xi_prop, xi_prop = xi_cur)
  
  # CURRENT 
  # Determine q_beta
  par_c1 = betapar(tsq = tau[m_star], sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_cur[m_star]+1):(xi_cur[m_star+1])], B = B)
  par_c2 = betapar(tsq = tau[m_star+1], sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_cur[m_star+1]+1):(xi_cur[m_star+2])], B = B)
  
  #extract beta_c1 and beta_c2
  beta_c1 = Beta[m_star, ]
  beta_c2 = Beta[m_star+1, ]
  
  # get q_betas
  q_betac1 = mvtnorm::dmvnorm(beta_c1, mean = par_c1$mean, sigma = par_c1$cov,log = TRUE)
  q_betac2 = mvtnorm::dmvnorm(beta_c2, mean = par_c2$mean, sigma = par_c2$cov, log = TRUE)
  
  # PROPOSAL
  # Determine q_beta
  par_p1 = betapar(tsq = tau[m_star], sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_prop[m_star]+1):(xi_prop[m_star+1])], B = B)
  par_p2 = betapar(tsq = tau[m_star+1], sigmasalpha = sigmasalpha, D = D, ts_seg = timeseries[(xi_prop[m_star+1]+1):(xi_prop[m_star+2])], B = B)
  
  # get proposal betas
  beta_p1 = c(mvtnorm::rmvnorm(1, mean = par_p1$mean, sigma = par_p1$cov))
  beta_p2 = c(mvtnorm::rmvnorm(1, mean = par_p2$mean, sigma = par_p2$cov))
  
  # get q_betas
  q_betap1 = mvtnorm::dmvnorm(beta_p1, mean = par_p1$mean, sigma = par_p1$cov, log = TRUE)
  q_betap2 = mvtnorm::dmvnorm(beta_p2, mean = par_p2$mean, sigma = par_p2$cov, log = TRUE)
  
  # TAU
  # proposal
  tau_p1 = 1/rgamma(1, (B + nu_0)/2, sum((beta_p1 * beta_p1)[-1] / D) / 2 + nu_0)
  tau_p2 = 1/rgamma(1, (B + nu_0)/2, sum((beta_p2 * beta_p2)[-1] / D) / 2 + nu_0)
  
  # current
  tau_c1 = tau[m_star]
  tau_c2 = tau[m_star + 1]
  
  # Likelihood
  # proposal:
  L_p1 = log_likelihood_adapt(b=beta_p1, Psi = par_p1$Psi, sumPsi = par_p1$sumPsi, perio = par_p1$perio) 
  L_p2 = log_likelihood_adapt(b=beta_p2, Psi = par_p2$Psi, sumPsi = par_p2$sumPsi, perio = par_p2$perio)
  # current
  L_c1 = log_likelihood_adapt(b=beta_c1, Psi = par_c1$Psi, sumPsi = par_c1$sumPsi, perio = par_c1$perio) 
  L_c2 = log_likelihood_adapt(b=beta_c2, Psi = par_c2$Psi, sumPsi = par_c2$sumPsi, perio = par_c2$perio)
  
  # Posterior
  # proposal
  post_prop = L_p1 + L_p2 +
    dxi(xi = xi_prop, tmin = tmin) + 
    sum(dnorm(beta_p1,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_p1)), log = TRUE)) +
    sum(dnorm(beta_p2,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_p2)), log = TRUE)) +
    invgamma::dinvgamma(tau_p1, nu_0/2, nu_0, log = TRUE) +
    invgamma::dinvgamma(tau_p2, nu_0/2, nu_0, log = TRUE)
  
  # current
  post_cur = L_c1 + L_c2 +
    dxi(xi = xi_cur, tmin = tmin) +
    sum(dnorm(beta_c1,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_c1)), log = TRUE)) +
    sum(dnorm(beta_c2,mean = 0, sd = sqrt(c(sigmasalpha, D * tau_c2)), log = TRUE)) +
    invgamma::dinvgamma(tau_c1, nu_0/2, nu_0, log = TRUE) +
    invgamma::dinvgamma(tau_c2, nu_0/2, nu_0, log = TRUE)
  
  # Calculate acceptance ratio
  q_prop = log(q_xip) + q_betap1 + q_betap2 + 
    invgamma::dinvgamma(tau_p1, (B + nu_0)/2, sum((beta_p1 * beta_p1)[-1] / D) / 2 + nu_0, log = TRUE)+
    invgamma::dinvgamma(tau_p2, (B + nu_0)/2, sum((beta_p2 * beta_p2)[-1] / D) / 2 + nu_0, log = TRUE)
  
  q_cur = log(q_xic) + q_betac1 + q_betac2 + 
    invgamma::dinvgamma(tau_c1, (B + nu_0)/2, sum((beta_c1 * beta_c1)[-1] / D) / 2 + nu_0, log = TRUE)+
    invgamma::dinvgamma(tau_c2, (B + nu_0)/2, sum((beta_c2 * beta_c2)[-1] / D) / 2 + nu_0, log = TRUE)
  
  
  A = exp(post_prop + q_cur - (post_cur + q_prop))
  # Metropolis Step
  U = runif(1)
  if(U < A){
    # Accept
    Beta_prop = Beta
    Beta_prop[c(m_star,m_star+1), ] = rbind(beta_p1,beta_p2)
    tau_prop = tau
    tau_prop[c(m_star,m_star+1)] = c(tau_p1,tau_p2)
    return(list("Beta" = Beta_prop, "xi" = xi_prop, "tau" = tau_prop))
  }else{
    # Reject
    return(list("Beta" = Beta, "xi" = xi_cur, "tau" = tau))
  }
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
}





