#' Whittle Post
#'
#' @param ab  : (K+1) column vector of alpha and beta values
#' @param X  : basis vector matrix
#' @param sumX : column sums of the X matrix
#' @param tsq   : tau squared scalar
#' @param perio : log of the periodogram
#' @param sigmasalpha : prior intercept variance, variance associated with the alpha_0 prior
#' @param D : vector of length K with strictly positive entries
#'
#' @return
#' @export
#'
#' @examples
whittle_post <- function(ab, X,sumX, tsq, perio, sigmasalpha, D){
  # pull out alpha
  a = ab[1]
  # pull out beta
  b = ab[-1]
  -1/2 * (sum(b * b / D) / tsq + a^2 / sigmasalpha + nrow(X) * a + crossprod(sumX, b) +
            sum(exp(perio - a - X %*% b)))
}