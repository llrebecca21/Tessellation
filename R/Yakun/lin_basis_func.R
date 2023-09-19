# Variables
# freq : ((nfreq + 1) x 1) vector of frequency values
# nbeta : numeric singular value
lin_basis_func <- function(freq, nbeta)
  {
    n <- length(freq)
    omega <- matrix(1, n, nbeta)
    for (j in 2:nbeta)
    {
      omega[, j] <- sqrt(2) * cos((j-1) * pi * freq)/(pi*(j-1))
    }
    
    return(omega)
  }
