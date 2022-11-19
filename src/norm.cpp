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


//[[Rcpp::export]]
double norm(arma::vec x)
{
  
  arma::vec m=x%x;
  
  return sqrt(sum(m));
}
