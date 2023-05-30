#' Beta Conditional Posterior for mulitple time series
#'
#' @param b : beta
#' @param X : basis function matrix
#' @param sumX : column sums of the X matrix
#' @param tsq  : tau squared scalar
#' @param y_bar : average periodogram vector
#' @param D : vector of length K with strictly positive entries
#'
#' @return
#' @export
#'
#' @examples
posterior_multiple <- function(b, X, sumX, tsq, y_bar, D, num_timeseries){
  -num_timeseries/2 * (sum(b * b / D) / (num_timeseries) + crossprod(sumX, b) +
            sum(y_bar / exp(X %*% b)))
}