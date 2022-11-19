#include <cmath>
#include <cstddef>
#include <algorithm>

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <boost/math/tools/roots.hpp>


//[[Rcpp::depends(BH)]]
//[[Rcpp::depends(RcppArmadillo)]]

namespace tools = boost::math::tools;
using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::uvec smallmove_distanceC(arma::mat x, arma::vec index, int Sindex, arma::vec w)
{
  int ii,zz;
  int n=index.n_elem;
  
  
  
  
  arma::vec dists=arma::zeros<arma::vec>(n);
  arma::uvec ind=arma::zeros<arma::uvec>(n);
  
  
  for(ii=0;ii<n;ii++)
  {
    
    arma::vec distance=arma::zeros<arma::vec>(x.n_cols);
    
    
    for(zz=0;zz<x.n_cols;zz++)
    {
      distance(zz)=pow((x(index(ii)-1,zz)-x(Sindex-1,zz)),2);
    }
    
    
    dists(ii)=sum(pow(w,2)%distance);
    
    
  }
  
  ind=sort_index(dists)+1;
  
  
  //return(distance);
  return(ind);
}
