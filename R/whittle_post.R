#' Title
#'
#' @param ab  : (K+1) column vector of alpha and beta values
#' @param sumX : row sums of the X matrix
#' @param t   : scalar containing tau value
#' @param perio : log of the periodogram
#' @param sigmasalpha : prior intercept variance, variance associated with the alpha_0 prior
#'
#' @return
#' @export
#'
#' @examples
whittle_post <- function(ab, sumX, t, perio, sigmasalpha){
  # pull out alpha
  a = ab[1]
  # pull out beta
  b = ab[-1]
  -1/2 * (crossprod(b) / t + a^2 / sigmasalpha + nrow(X) * a + sumX %*% b +
            sum(exp(perio - a - X %*% b)))
}