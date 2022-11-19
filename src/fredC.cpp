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
double fredC(double beep, arma::vec gq, arma::vec trust_beta,double C1,double C2,double r)
{
  double result;
  
  if(beep==0)
  {
    if(C2>0)
    {
      result= -1.0/r;
      
    }else{
      
      result=sqrt(1.0 / C1) - 1.0 / r;
    }
    
  }else{
    
    result=sqrt(1.0 / sum(pow((gq / (trust_beta + beep)), 2))) - 1.0 / r;
  }
  
  return result;
}
