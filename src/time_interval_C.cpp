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
Rcpp::List time_interval_C(arma::mat x, arma::vec S, arma::vec prt, int NumObs)
{
  
  arma::mat interval=arma::zeros<arma::mat>(S.n_elem, NumObs);
  arma::mat tmin=arma::zeros<arma::mat>(S.n_elem, NumObs);
  arma::mat tmax=arma::zeros<arma::mat>(S.n_elem, NumObs);
  
  for(int i=0; i<S.n_elem; i++)
  {
    arma::uvec index=find((prt-1)==i);
    arma::mat x_s=x.rows(index);
    
    
    
    for(int j=0; j< NumObs; j++)
    {
      
      if(any((x_s.col(0)-1)==j)) // index is j
      {
        
        arma::vec temp_index=x_s.col(0);
        arma::uvec idx=find((temp_index-1)==j);
        arma::mat temp=x_s.rows(idx);
        arma::vec time=temp.col(1);
        
        
        tmin(i,j)=min(time);
        tmax(i,j)=max(time);
        interval(i,j)=tmax(i,j)-tmin(i,j)+1;
        
        
      }else{
        
        
        tmin(i,j)=-1;
        tmax(i,j)=-1;
        interval(i,j)=-1;
        
      }
    }
    
  }
  
  Rcpp::List re=List::create(Named("interval",0), _["tmin"]=0, _["tmax"]=0);
  
  re[0]=interval;
  re[1]=tmin;
  re[2]=tmax;
  
  
  return(re);
  
  
}
