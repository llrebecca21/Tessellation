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

// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat const &x, 
                       arma::vec const &mean, 
                       arma::mat const &sigma, 
                       bool const logd = false) { 
  arma::vec const distval = Mahalanobis(x,  mean, sigma);
  double const logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec const logretval = 
    -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
    
    if (logd)
      return logretval;
    return exp(logretval);
}
