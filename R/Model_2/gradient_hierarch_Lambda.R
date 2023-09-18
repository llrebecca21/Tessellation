#' gradient_hierarch_Lambda
#'
#' @param br 
#' @param b 
#' @param sumPsi 
#' @param Psi 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
gradient_hierarch_Lambda = function(br, b, sumPsi, Psi, y, lambda){
  -lambda %*% (br - b) - sumPsi + colSums(Psi * c(y / exp(Psi %*% br)))
}