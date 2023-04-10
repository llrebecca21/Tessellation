source("R/lin_basis_func.R")

#' gr
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
gr <- function(param,
               index,
               nseg_time_temp,
               y,
               tau_temp,
               nbeta,
               nbasis,
               sigmasqalpha) {
  g <- matrix(0, nbeta, 1)
  g[1, ] <- param[1] / sigmasqalpha
  g[2:nbeta, ] <- param[2:nbeta] / tau_temp



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

      temp_mat <- nu_mat[2:(n1 + 1), ] * repmat((1 - exp(as.matrix(ydev[2:(n1 + 1)]))), 1, dim(nu_mat)[2])
      g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1])) %*%
        nu_mat[1, ])
    } else {
      nfreq <- floor(nseg_time_temp[l] / 2)
      freq <- (0:nfreq) / (2 * nfreq)
      nu_mat <- lin_basis_func(freq, nbeta)
      ydev <- y[[uu]] - nu_mat %*% param ### ydev is a vector

      n1 <- nfreq

      temp_mat <- nu_mat[2:n1, ] * repmat((1 - exp(as.matrix(ydev[2:n1]))), 1, dim(nu_mat)[2])
      g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1])) %*% nu_mat[1, ]) +
        t(0.5 * (1 - exp(ydev[n1 + 1])) %*% nu_mat[n1 + 1, ])
    }
  }

  return(g)
}
