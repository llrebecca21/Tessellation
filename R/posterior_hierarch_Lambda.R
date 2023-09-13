#' posterior_hierarch_Lambda : 
#'
#' @param br 
#' @param b 
#' @param Psi 
#' @param sumPsi 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
posterior_hierarch_Lambda <- function(br, b, Psi, sumPsi, y){
  -(1/2) * crossprod(br - b, lambda %*% (br - b))  - crossprod(sumPsi, br) -
    sum(y / exp(Psi %*% br))
}