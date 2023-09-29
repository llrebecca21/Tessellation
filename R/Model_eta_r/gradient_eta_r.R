#' function that calculates the gradient of the conditional posterior distribution of the beta^r
#'
#' @param br 
#' @param sumPsi 
#' @param Psi 
#' @param y 
#' @param Sigma 
#' @param eta_r 
#'
#' @return
#' @export
#'
#' @examples
gradient_eta_r = function(br, sumPsi, Psi, y, Sigma, eta_r){
  -sumPsi + colSums(Psi * c(y / exp(Psi %*% br))) - 2 * br * eta_r / Sigma
}