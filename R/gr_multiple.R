#' Gradient Function for multiple stationary time series with same underlying spectrum
#'
#' @param b 
#' @param X 
#' @param sumX 
#' @param tsq 
#' @param y_bar 
#' @param D 
#' @param num_timeseries 
#'
#' @return
#' @export
#'
#' @examples
gr_multiple <- function(b, X, sumX, tsq, y_bar, D, num_timeseries){
  # derivative of log of the whittle_post with respect to beta (vector)
  derivbeta = -(num_timeseries/2) * (2 * (b / D) / (num_timeseries) + sumX - colSums(X * c(y_bar / exp(X %*% b))))
}