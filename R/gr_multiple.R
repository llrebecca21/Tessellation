#' Gradient Function for multiple stationary time series with same underlying spectrum
#'
#' @param b 
#' @param X 
#' @param sumX 
#' @param tsq 
#' @param y_bar 
#' @param D 
#' @param B 
#'
#' @return
#' @export
#'
#' @examples
gr_multiple <- function(b, Psi, sumPsi, tsq, y_bar, D, R){
  # derivative of log of the whittle_post with respect to beta (vector)
  derivbeta = -R * ((b / D) / R + sumPsi - colSums(Psi * c(y_bar / exp(Psi %*% b))))
}