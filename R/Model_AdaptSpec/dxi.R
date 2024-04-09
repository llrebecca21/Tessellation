# Function that calculates the distribution of the xi (xi's prior)

dxi = function(xi,tmin){
  S = length(xi)-1
  m = 2:S
  # return log of dxi
  return(-1*sum(log((xi[S+1] - xi[m-1] - (S - m + 1)*tmin + 1))))
}


