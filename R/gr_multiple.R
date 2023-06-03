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
gr_multiple <- function(b, X, sumX, tsq, y_bar, D, B){
  # derivative of log of the whittle_post with respect to beta (vector)
  derivbeta = -(B/2) * (2 * (b / D) / (B) + sumX - colSums(X * c(y_bar / exp(X %*% b))))
}