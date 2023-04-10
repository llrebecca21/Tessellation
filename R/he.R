source("R/lin_basis_func.R")
#' he
#'
#' @param param :
#' @param index :
#' @param nseg_time_temp : 
#' @param y :
#' @param tau_temp : 
#' @param nbeta :
#' @param nbasis :
#' @param sigmasqalpha :
#'
#' @return
#' @export
#'
#' @examples
he <- function(param,
               index,
               nseg_time_temp,
               y,
               tau_temp,
               nbeta,
               nbasis,
               sigmasqalpha) {
  h <- matrix(0, nbeta, nbeta)
  h[1, 1] <- 1 / sigmasqalpha
  h[2:nbeta, 2:nbeta] <- 1 / tau_temp * diag(nbasis)


  uu <- 0
  for (l in index)
  {
    uu <- uu + 1
    if (nseg_time_temp[l] %% 2 == 1) {
      nfreq <- floor(nseg_time_temp[l] / 2)
      freq <- (0:nfreq) / (2 * nfreq)
      nu_mat <- lin_basis_func(freq, nbeta)
      ydev <- y[[uu]] - nu_mat %*% param ### ydev is a vector

      n1 <- nfreq

      jj <- seq(1:nbeta)
      big_mat <- repmat(t(nu_mat[2:(n1 + 1), ]), nbeta, 1) * t(nu_mat[2:(n1 + 1), repmat(jj, nbeta, 1)])

      coef_mat <- repmat(exp(t(ydev[2:(n1 + 1)])), nbeta^2, 1)
      h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1)), c(1, 2), sum) +
        t(0.5 * exp(ydev[1]) %*% nu_mat[1, ]) %*% nu_mat[1, ]
    } else {
      nfreq <- floor(nseg_time_temp[l] / 2)
      freq <- (0:nfreq) / (2 * nfreq)
      nu_mat <- lin_basis_func(freq, nbeta)
      ydev <- y[[uu]] - nu_mat %*% param ### ydev is a vector

      n1 <- nfreq

      jj <- seq(1:nbeta)
      big_mat <- repmat(t(nu_mat[2:n1, ]), nbeta, 1) * t(nu_mat[2:n1, repmat(jj, nbeta, 1)])

      coef_mat <- repmat(exp(t(ydev[2:n1])), nbeta^2, 1)
      h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1 - 1)), c(1, 2), sum) +
        t(0.5 * exp(ydev[1]) %*% nu_mat[1, ]) %*% nu_mat[1, ] +
        t(0.5 * exp(ydev[n1 + 1]) %*% nu_mat[n1 + 1, ]) %*% nu_mat[n1 + 1, ]
    }
  }

  return(h)
}
