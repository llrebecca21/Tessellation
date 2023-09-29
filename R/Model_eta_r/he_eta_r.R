#' function that calculates the hessian of the conditional posterior distribution of beta^r
#'
#' @param Psi 
#' @param y 
#' @param br 
#' @param eta_r 
#' @param Sigma 
#'
#' @return
#' @export
#'
#' @examples
he_eta_r = function(Psi, y, br, eta_r, Sigma){
  he1 = -crossprod(Psi , Psi * c(y / exp(Psi %*% br)))
  diag(he1) = diag(he1) - 2 * eta_r / Sigma
  return(he1)
}