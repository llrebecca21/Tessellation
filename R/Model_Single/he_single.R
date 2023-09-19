#' Calculate the negative hessian for a single stationary time series
#'
#' @param ab : (K+1) column vector of alpha and beta values
#' @param X  : basis vector matrix
#' @param sumX : column sums of the X matrix
#' @param tsq  : tau squared scalar
#' @param perio : log of the periodogram
#' @param sigmasalpha : prior intercept variance, variance associated with the alpha_0 prior
#' @param D : vector of length K with strictly positive entries
#'
#' @return the negative of the hessian
#' @export
#'
#' @examples
he_single <- function(ab, X, sumX, tsq, perio, sigmasalpha, D){
  # pull out alpha
  a = ab[1]
  # pull out beta
  b = ab[-1]
  # Calculate the residual vector (equation inside the exp)
  exp_res = c(exp(perio - X %*% b - a))
  # he_aa : stands for second partial derivative both with respect to alpha
  he_aa <- 1/sigmasalpha + 0.5 * sum(exp_res)
  # he_ab : stands for partial with respect to alpha and then beta
  Diag_res_X = X * exp_res
  he_ab = 0.5 * colSums(Diag_res_X)
  # he_bb : stands for second partial derivative both with respect to beta (vector)
  he_bb = 0.5 * crossprod(X, Diag_res_X)
  # update the diagonal values of he_bb only
  diag(he_bb) = diag(he_bb) + 1 / (tsq * D)
  he_matrix = matrix(NA, nrow = length(ab), ncol = length(ab))
  he_matrix[1,1] <- he_aa
  he_matrix[1,-1] <- t(he_ab)
  he_matrix[-1,1] <- he_ab
  he_matrix[-1,-1] <- he_bb
  return(he_matrix)
}