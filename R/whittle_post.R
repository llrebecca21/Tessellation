#' Title
#'
#' @param ab   : (K + 1 x 1) vector that stores alpha and beta values
#' @param sumX : Squared Euclidean norm of the matrix X
#' @param t    : vector that stores tau values
#'
#' @return
#' @export
#'
#' @examples
whittle_post <- function(ab, sumX, t){
  # pull out alpha
  a = ab[1]
  # pull out beta
  b = ab[-1]
  -1/2 * (crossprod(b) / t + a^2 / sigmasalpha + nrow(X) * a + sumX %*% b +
            sum(exp(perio - a - X %*% b)))
}