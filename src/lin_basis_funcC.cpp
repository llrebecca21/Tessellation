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
arma::mat lin_basis_funcC(arma::rowvec freq, int nbeta)
{
  
  int n=freq.n_elem;
  arma::mat out;
  out.ones(n,nbeta);
  
  int j;
  
  for(j=1; j<nbeta; j++)
  {
    
    out.col(j)=sqrt(2) * cos(j * datum::pi * freq.t())/(datum::pi*j);
    
  }
  
  
  return(out);
  
}
