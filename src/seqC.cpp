#include <cmath>
#include <cstddef>
#include <algorithm>

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <omp.h>


#include <boost/math/tools/roots.hpp>


//[[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::plugins(openmp)]]


namespace tools = boost::math::tools;
using namespace Rcpp;
using namespace arma;

// matrix_big_mat internal
// purpose of this function: 

//[[Rcpp::export]]
arma::rowvec seqC(double start, double end, double nfreq)
{
  arma::rowvec re(nfreq+1,fill::zeros);
  
  
  int i;
  for(i=1; i<(nfreq+1); i++)
  {
    re(i)=re(i-1)+1/(nfreq*2);
  }
  
  return(re);
  
}
