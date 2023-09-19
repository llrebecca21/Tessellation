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
double whittle_like_C(Rcpp::List ytemp, arma::vec index, arma::vec nseg_time_temp, arma::vec beta, int nbeta, int nbasis)
{
  
  double log_prop_spec_dens=0;
  
  
  for (int uu=0;uu<index.n_elem; uu++)
  {
    
    
    int l=index(uu)-1;
    
    
    int nfreq;
    nfreq=std::floor(nseg_time_temp(l)/2);
    
    
    arma::rowvec freq=seqC(0.0,0.5,nfreq);
    arma::mat nu_mat=lin_basis_funcC(freq,nbeta);
    
    
    arma::mat fhat=nu_mat*beta;  // nfreq+1 by 1
    arma::vec temp_y=ytemp[uu];
    
    
    
    int n1=nfreq;
    int n=nseg_time_temp(l);
    int mod=n-n1*2;
    
    
    
    if (mod==1)
    {
      
      log_prop_spec_dens=log_prop_spec_dens-sum(fhat.rows(1,nfreq)+exp(temp_y.subvec(1,nfreq)-fhat.rows(1,nfreq)))-
        0.5*(fhat(0)+exp(temp_y(0)-fhat(0)))-0.5*nfreq*log(2*datum::pi);
      
      
    }else
    {
      
      log_prop_spec_dens=log_prop_spec_dens-sum(fhat.rows(1,nfreq-1)+exp(temp_y.subvec(1,nfreq-1)-fhat.rows(1,nfreq-1)))-
        0.5*(fhat(0)+exp(temp_y(0)-fhat(0)))-0.5*(fhat(nfreq)+exp(temp_y(nfreq)-fhat(nfreq)))-
        0.5*nfreq*log(2*datum::pi);
      
      
    }
    
  }
  
  
  
  
  return(log_prop_spec_dens);
  
}
