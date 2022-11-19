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
Rcpp::List postbeta_C(arma::vec index,arma::vec tmax, arma::vec tmin, arma::vec nseg_time_temp, arma::mat x_t, double  tau_temp, int nbeta, int nbasis, double sigmasqalpha)
{
  
  Rcpp::List y(index.n_elem);
  Rcpp::List yy(index.n_elem);
  
  for (int uu=0; uu<index.n_elem; uu++)
  {
    
    
    int l=index(uu)-1;
    
    int nfreq;
    nfreq=std::floor(nseg_time_temp(l)/2);
    
    y[uu]=log(pow(abs(fft(x_t(span(tmin(l)-1,tmax(l)-1), l))) , 2) / nseg_time_temp(l));
    
    
    arma::vec ytemp=y[uu];
    yy[uu]=ytemp.rows(0,nfreq);
    
    
  }
  
  
  arma::vec param(nbeta,fill::zeros);
  double rinit=1;
  double rmax=100;
  arma::vec parscale(nbeta,fill::ones);
  int iterlim=100;
  double fterm=0.0000000149;
  double mterm=0.0000000149;
  bool minimize=true;
  
  
  Rcpp::List post=trust_C(param, rinit, rmax, parscale, iterlim, fterm, mterm, minimize, index, nseg_time_temp, yy, tau_temp, nbeta, nbasis, sigmasqalpha);
  arma::mat var=post["hessian"];
  
  
  Rcpp::List re=List::create(Named("beta_mean",0), _["beta_var"]=0, _["y"]=0);
  
  re[0]=post["argument"];
  re[1]=inv(var);
  re[2]=y;
  
  
  
  return(re);
  
  
}  
