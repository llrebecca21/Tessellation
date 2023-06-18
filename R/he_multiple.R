#' Calculate the negative Hessian of beta for multiple stationary time series with same underlying spectrum
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
he_multiple <- function(b, Psi, y_bar, D, R){
  he_bb = R * crossprod(Psi , Psi * c(y_bar / exp(Psi %*% b)))
  diag(he_bb) =  1 / (D) + diag(he_bb)
  return(he_bb)
}