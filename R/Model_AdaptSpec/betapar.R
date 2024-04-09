# create a function to get Gaussian parameters used to estimate a beta proposal
betapar = function(tsq, sigmasalpha, D, ts_seg, B){
  n = length(ts_seg)
  J = floor((n - 1) / 2)
  omega = (2 * pi * (0:J)) / n
  Psi = cbind(1, outer(X = omega, Y = 1:B, FUN = function(x,y){cos(y * x) * sqrt(2)}))
  sumPsi = c(crossprod(rep(1, nrow(Psi)), Psi))
  perio = log((abs(fft(ts_seg))^2 / n))
  # Maximum A Posteriori (MAP) estimate : finds the alpha and beta that gives us the mode of the conditional posterior of beta and alpha_0 conditioned on y
  map = optim(par = b, fn = beta_cond_post, gr = gr_adapt, method ="BFGS", control = list(fnscale = -1),
              Psi = Psi, sumPsi = sumPsi, tsq = tsq, perio = perio, Sigma = c(sigmasalpha, D * tsq))$par
  # Call the hessian function
  norm_precision = -1 * he_adapt(beta_mstar = map, Psi = Psi, sumPsi = sumPsi, Sigma =  c(sigmasalpha, D * tsq), perio = perio)
  return(list("mean" = map, "cov" = solve(norm_precision), "perio" = perio, "Psi" = Psi, "sumPsi" = sumPsi))
}