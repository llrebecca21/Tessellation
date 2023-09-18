#' he_hierarch_Lambda
#'
#' @param br 
#' @param lambda 
#' @param Psi 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
he_hierarch_Lambda = function(br, lambda, Psi, y){
  -lambda -  crossprod(Psi , Psi * c(y / exp(Psi %*% br)))
}