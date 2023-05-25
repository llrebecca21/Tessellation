#' Calculate the negative Hessian of beta for multiple stationary time series with same underlying spectrum
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
he_multiple <- function(b, X, sumX, tsq, y_bar, D, num_timeseries){
  he_bb = (num_timeseries / 2) * crossprod(X , X * c(y_bar / exp(X %*% b)))
  diag(he_bb) =  1 / (tsq * D) + diag(he_bb)
  return(he_bb)
}