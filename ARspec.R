#' ARspec: Calculates the true spectral density for an AR process
#'
#' @param phi 
#' @param sigma 
#' @param freq 
#'
#' @return
#' @export
#'
#' @examples
ARspec <- function(phi, sigma, freq) {
  dim <- dim(phi)
  len <- dim[2]
  phispect <- array(0, c(2, 2, length(freq)))
  spect <- array(0, c(2, 2, length(freq)))

  for (k in 1:length(freq))
  {
    phispect[, , k] <- diag(2)
    for (j in 1:(len / 2))
    {
      if (j == 1) {
        bigmat <- phi[, 1:2] * exp(-2 * pi * sqrt(-1 + 0i) * freq[k])
      } else {
        bigmat <- phi[, (2 * j - 1):(2 * j)] * exp(-2 * j * pi * sqrt(-1 + 0i) * freq[k])
      }

      phispect[, , k] <- phispect[, , k] - bigmat
    }
    spect[, , k] <- sigma %*% solve(phispect[, , k]) %*% Conj(solve(phispect[, , k]))
  }

  return(spect)
}
