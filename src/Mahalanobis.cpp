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
arma::vec Mahalanobis(arma::mat const &x, 
                      arma::vec const &center, 
                      arma::mat const &cov) {
  arma::mat x_cen = x.t();
  x_cen.each_col() -= center;
  arma::solve(x_cen, arma::trimatl(chol(cov).t()), x_cen);
  x_cen.for_each( [](arma::mat::elem_type& val) { val = val * val; } );
  return arma::sum(x_cen, 0).t();    
}
