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
//[[Rcpp::export]]
arma::mat matrix_big_mat(arma::mat a, arma::mat b)
{
  arma::mat out;
  out.zeros(a.n_rows,b.n_rows*b.n_cols);
  
  int i,j,k;
  
  
  for(i=0;i<a.n_rows;i++)
  {
    for(j=0;j<a.n_cols;j++)
    {
      for(k=0;k<a.n_cols;k++)
      {
        
        out(i,(j*a.n_cols+k))=a(i,j);
        
      }
    }
  }
  
  return(out);
}

