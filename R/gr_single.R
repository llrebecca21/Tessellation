#' Gradient Function for a single stationary time series
#'
#' @param ab : (K+1) column vector of alpha and beta values
#' @param X  : basis vector matrix
#' @param sumX : column sums of the X matrix
#' @param tsq  : tau squared scalar
#' @param perio : log of the periodogram
#' @param sigmasalpha : prior intercept variance, variance associated with the alpha_0 prior
#' @param D : vector of length K with strictly positive entries
#'
#' @return
#' @export
#'
#' @examples
gr_single <- function(ab, X, sumX, tsq, perio, sigmasalpha, D){
  # pull out alpha
  a = ab[1]
  # pull out beta
  b = ab[-1]
  # Calculate the residual vector (equation inside the exp)
  exp_res = c(exp(perio - X %*% b - a))
  # derivative of log of the whittle_post with respect to alpha
  derivalpha = -0.5 * (2 * a / sigmasalpha - length(sumX) - sum(exp_res))
  # derivative of log of the whittle_post with respect to beta (vector)
  derivbeta = -0.5 * (2 * (b / D) / tsq - sumX - colSums(X * exp_res))
  # combine the two derivative for the gradient vector
  grad = c(derivalpha, derivbeta)
  return(grad)
}