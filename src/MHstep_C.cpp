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
Rcpp::List MHstep_C(arma::mat x, arma::mat x_t, Rcpp::List TT, int Ntime, int Mmax, int Tmin, int NumObs, double var_w, int nbeta, int nbasis, double sigmasqalpha, int ncores)
{
  
  
  
  int TTM=TT["M"];
  arma::vec accept(TTM,fill::zeros);
  arma::vec alpha(TTM,fill::zeros);
  
  
  // proposed partition
  arma::mat x_temp=x.cols(2,x.n_cols-1);
  arma::vec prt=distance_partitionC(x_temp,TT["S"],TT["w"]);
  
  // time interval
  Rcpp::List interval_curr=time_interval_C(x, TT["S"], prt, NumObs);
  arma::mat tmax=interval_curr["tmax"];
  arma::mat tmin=interval_curr["tmin"];
  arma::mat interval=interval_curr["interval"];
  
  
  // beta, tau, and g
  arma::vec tau_prop=TT["tau"];
  arma::vec g_prop=TT["g"];
  arma::mat beta=TT["beta"];
  
  
  // self defined basis
  int nfreq=50;
  arma::rowvec freq=seqC(0.0,0.5,nfreq);
  arma::mat nu_mat=lin_basis_funcC(freq,nbeta);
  
  
  
#pragma omp parallel for num_threads(1)
  
  for(int i=0; i<TTM; i++)
  {
    
    arma::uvec index=find((prt-1)==i);
    arma::mat prop_obs=x.rows(index);
    arma::vec prop_obs_index=unique(prop_obs.col(0)); // unique index column
    
    arma::vec curr_tmax=conv_to<vec>::from(tmax.row(i));
    arma::vec curr_tmin=conv_to<vec>::from(tmin.row(i));
    arma::vec curr_interval=conv_to<vec>::from(interval.row(i));
    double curr_tau=tau_prop(i);
    
    
    //optimization for beta
    Rcpp::List postbeta_max=postbeta_C(prop_obs_index, curr_tmax, curr_tmin, curr_interval, x_t, curr_tau, nbeta, nbasis, sigmasqalpha);
    arma::mat beta_var=postbeta_max["beta_var"];
    arma::vec curr_beta=conv_to<vec>::from(beta.row(i));
    
    // use current beta to generate new beta
    //arma::vec beta_prop=mvnrnd(curr_beta,0.5*(beta_var+beta_var.t()));
    
    // use max mean to generate new beta
    arma::vec beta_mean=postbeta_max["beta_mean"];
    arma::vec beta_prop=mvnrnd(beta_mean,0.5*(beta_var+beta_var.t()));
    
    // proposal probability for beta
    arma::mat beta_mat(beta_prop);
    //arma::vec log_beta_prop=dmvnorm_arma(beta_mat.t(),curr_beta,0.5*(beta_var+beta_var.t()),true);
    arma::vec log_beta_prop=dmvnorm_arma(beta_mat.t(),beta_mean,0.5*(beta_var+beta_var.t()),true);
    
    // prior probability for beta
    arma::mat zero_mean=arma::zeros<arma::mat>(nbeta, 1);
    arma::vec diag(nbeta, fill::ones);
    diag(0)=sigmasqalpha;
    arma::mat var=diagmat(diag);
    arma::vec log_beta_prior_prop=dmvnorm_arma(beta_mat.t(), zero_mean, var, true);
    
    // fhat and Whittle likelihood
    arma::mat fhat_prop=nu_mat*beta_mat;
    double log_prop_spec_dens=whittle_like_C(postbeta_max["y"],prop_obs_index,curr_interval,beta_prop,nbeta,nbasis);
    double loglike_prop=log_prop_spec_dens;
    
    
    // for alpha
    double log_proposal_prop = log_beta_prop(0);
    arma::vec log_beta=TT["log_beta"];
    double log_proposal_curr = log_beta(i);
    
    
    double beta_lprior_prop=log_beta_prior_prop(0);
    arma::vec beta_lprior=TT["beta_lprior"];
    double beta_lprior_curr=beta_lprior(i);
    
    
    double lprior_prop=beta_lprior_prop;
    double lprior_curr=beta_lprior_curr;
    
    
    
    // log-MH ratio
    arma::vec loglike=TT["loglike"];
    arma::vec ratio={1,exp(loglike_prop+lprior_prop-loglike(i)-lprior_curr+log_proposal_curr-log_proposal_prop)};
    alpha(i)=min(ratio);
    
    
    // create new TT
    // Rcpp::List TT_new = List::create(Named("S",TT["S"]), _["M"]=TT["M"], _["w"]=TT["w"], _["loglike"]=TT["loglike"],
    //                                  _["beta"]=TT["beta"], _["tau"]=TT["tau"],_["g"]=TT["g"],_["log_beta"]=TT["log_beta"],
    //                                    _["M_lprior"]=TT["M_lprior"],_["beta_lprior"]=TT["beta_lprior"],
    //                                   _["fhat_prop"]=TT["fhat_prop"],_["prt"]=TT["prt"] );                                               );
    
    double rr = randu();
    if(alpha(i)>rr)
    {
      
      loglike(i)=loglike_prop;
      TT["loglike"]=loglike;
      
      log_beta(i)=log_beta_prop(0);
      TT["log_beta"]=log_beta;
      
      beta_lprior(i)=beta_lprior_prop;
      TT["beta_lprior"]=beta_lprior;
      
      beta.row(i)=conv_to<mat>::from(beta_prop).t();
      TT["beta"]=beta;
      
      arma::mat fhat_prop_mat=TT["fhat_prop"];
      fhat_prop_mat.col(i)=fhat_prop;
      TT["fhat_prop"]=fhat_prop_mat;
      
      accept(i)=1;
      
    }else{
      
      alpha(i)=-1;
      accept(i)=0;
    }
    
  }
  
  
  Rcpp::List re = List::create(Named("alpha",0), _["TT"]=0, _["accept"]=0);
  
  
  re[0]=alpha;
  re[1]=TT;
  re[2]=accept;
  
  return(re);
  
  
  
  
}
