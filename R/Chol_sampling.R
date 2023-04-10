#' Chol_sampling
#'
#' @param Lt : (nbeta x nbeta) : Upper triangular matrix from Cholesky decomposition
#' @param d  : (numeric : nbeta) : Length of beta
#' @param beta_c : (nbeta x 1) : Vector of current beta values
#'
#' @return
#' @export
#'
#' @examples
Chol_sampling <- function(Lt, d, beta_c) {
  u <- solve(Lt, rnorm(n = d))
  beta_prop <- beta_c + u

  return(beta_prop)
}
