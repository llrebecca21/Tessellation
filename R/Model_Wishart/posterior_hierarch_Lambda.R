#' posterior_hierarch_Lambda : function calculates the conditional posterior distribution of bbmath{B}
#'
#' @param br :
#' @param b :
#' @param Psi : (J x (B+1)) 
#' @param sumPsi : 
#' @param y : 
#'
#' @return
#' @export
#'
#' @examples
posterior_hierarch_Lambda <- function(br, b, Psi, sumPsi, y, lambda){
  -(1/2) * crossprod(br - b, lambda %*% (br - b))  - crossprod(sumPsi, br) -
    sum(y / exp(Psi %*% br))
}