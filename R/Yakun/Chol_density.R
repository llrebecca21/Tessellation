#' Chol_density
#'
#' @param x  :
#' @param Lt :
#' @param d  :
#' @param m  :
#'
#' @return
#' @export
#'
#' @examples
Chol_density <- function(x, Lt, d, m) {
  x_tilde <- Lt %*% (x - m)
  # term1 : Leading normalizing constant
  term1 <- (-d / 2) * log(2 * pi)
  # term2 : determinant of Q
  term2 <- sum(log(diag(L)))
  # term3: log of the exp term of the normal dist
  term3 <- (-1 / 2) * sum(x_tilde^2)
  log_beta_prop <- term1 + term2 + term3

  return(log_beta_prop)
}
