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
arma::vec distance_partitionC(arma::mat x, arma::vec S, arma::vec w)
{
  int ii,jj,zz;
  int n=x.n_rows;
  int M=S.n_elem;
  
  
  arma::vec prt=arma::zeros<arma::vec>(n);
  arma::vec dists=arma::zeros<arma::vec>(M);
  arma::uvec sort=arma::zeros<arma::uvec>(M);
  
  
  for(ii=0;ii<x.n_rows;ii++)
  {
    int ind;
    
    
    for(jj=0;jj<M;jj++)
    {
      
      arma::vec distance=arma::zeros<arma::vec>(x.n_cols);
      
      
      for(zz=0;zz<x.n_cols;zz++)
      {
        distance(zz)=pow((x(ii,zz)-x(S(jj)-1,zz)),2);
      }
      
      
      dists(jj)=sum(pow(w,2)%distance);
      
    }
    
    
    sort=sort_index(dists);
    ind=sort(0);
    prt(ii)=ind+1;
    
  }
  
  //return(distance);
  return(prt);
}
