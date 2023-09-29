#' Function that calculates the conditional posterior distribution of a single beta^r term
#'
#' @param sumPsi 
#' @param br 
#' @param Psi 
#' @param y 
#' @param eta_r 
#' @param Sigma 
#'
#' @return
#' @export
#'
#' @examples
posterior_eta_r = function(sumPsi, br, Psi, y, eta_r, Sigma){
  -crossprod(sumPsi, br) - sum(y / exp(Psi %*% br)) - sum(eta_r * br * br / Sigma)
}