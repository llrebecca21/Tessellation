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
Rcpp::List unirootC(double betadn, double betaup, arma::vec gq, arma::vec trust_beta,double C1,double C2, double r){
  
  Function f("uniroot"); 
  
  Function g("fredR");
  
  //return g(beep, _["gq"]=gq, _["beta"]=trust_beta);
  SEXP re=f(g, Named("lower")=betadn, _["upper"]=betaup,_["tol"]=0.000001, _["gq"]=gq, _["beta"]=trust_beta, _["C1"]=C1, _["C2"]=C2,_["r"]=r);
  
  //Rcpp::List result(re);
  return re;
  
}
