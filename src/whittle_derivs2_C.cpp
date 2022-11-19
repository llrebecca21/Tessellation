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
Rcpp::List whittle_derivs2_C(arma::vec x, arma::vec index, arma::vec nseg_time_temp, Rcpp::List ytemp, double tau_temp,int nbeta,int nbasis,double sigmasqalpha)
{
  
  
  arma::mat param(x); // nbeta by 1
  
  
  // f
  double f=0;
  
  
  // g
  arma::vec gr(nbeta,fill::zeros);
  
  gr(0)=param(0,0)/sigmasqalpha;
  gr.subvec(1,nbeta-1)=param.rows(1,nbeta-1)/tau_temp;
  arma::mat temp_mat;
  
  
  
  // h
  arma::mat he;
  he.zeros(nbeta, nbeta);
  
  he(0,0)=1.0/sigmasqalpha;
  he.submat(1,1,nbeta-1,nbeta-1)=(1.0/tau_temp)*(arma::eye(nbasis,nbasis));
  
  arma::mat jj=arma::linspace(1,nbeta,nbeta).t();
  arma::mat coef_mat;
  arma::mat big_mat;
  
  for (int uu=0;uu<index.n_elem; uu++)
  {
    
    
    int l=index(uu)-1;
    
    
    int nfreq;
    nfreq=std::floor(nseg_time_temp(l)/2);
    
    
    arma::rowvec freq=seqC(0.0,0.5,nfreq);;
    arma::mat nu_mat=lin_basis_funcC(freq,nbeta);
    
    
    arma::mat temp_numat=nu_mat*param;  // nfreq+1 by 1
    arma::vec temp_y=ytemp[uu];
    arma::vec ydev(nfreq+1,fill::zeros);
    
    // create ydev
    for(int i=0; i<(nfreq+1); i++)
    {
      
      ydev(i)=temp_y(i)-temp_numat(i);
      
    }
    
    
    
    
    int n1=nfreq;
    int n=nseg_time_temp(l);
    int mod=n-n1*2;
    
    
    
    if (mod==1)
    {
      // f
      arma::mat pro=(param.rows(1,nbeta-1).t()*param.rows(1,nbeta-1))/tau_temp+pow(param(0,0),2)/sigmasqalpha;
      arma::mat sym=nu_mat.row(0)*param;
      f=sum(nu_mat.rows(1,n1)*param + exp(ydev.subvec(1,n1)))+
        0.5*(sum(sym(0,0)+exp(ydev(0))))+
        0.5*pro(0,0);
      
      
      // g
      temp_mat=nu_mat.rows(1,n1)%repmat((1-exp(ydev.subvec(1,n1))),1,nu_mat.n_cols);
      gr=gr+sum(temp_mat,0).t()+(0.5*(1-exp(ydev(0)))*nu_mat.row(0)).t();
      
      
      // h
      big_mat=repmat(nu_mat.rows(1,n1).t(),nbeta,1) % (matrix_big_mat(nu_mat.rows(1,n1),repmat(jj,nbeta,1)).t());
      
      coef_mat=repmat(exp(ydev.subvec(1,n1).t()),pow(nbeta,2),1);
      
      arma::mat temp=big_mat % coef_mat;
      arma::cube c(nbeta, nbeta, n1);
      
      int jj=0;
      
      for(int z=0;z<c.n_slices;z++)
      {
        for(int j=0;j<c.n_cols;j++)
        {
          for(int q=0;q<c.n_rows;q++)
          {
            c(q,j,z)=temp(jj);
            jj=jj+1;
            
          }
        }
      }
      
      arma::mat sum_mat=sum(c,2);
      he=he+sum_mat+(0.5*exp(ydev(0))*nu_mat.row(0)).t() * nu_mat.row(0);
      
      
    }else
    {
      // f
      arma::mat pro=(param.rows(1,nbeta-1).t()*param.rows(1,nbeta-1))/tau_temp+pow(param(0,0),2)/sigmasqalpha;
      arma::mat sym=nu_mat.row(0)*param;
      arma::mat sym_one=nu_mat.row(n1)*param;
      f=sum(nu_mat.rows(1,n1-1)*param +exp(ydev.subvec(1,n1-1)))+
        0.5*(sum(sym(0,0)+exp(ydev(0))))+
        0.5*(sum(sym_one(0,0)+exp(ydev(n1))))+
        0.5*pro(0,0);
      
      // g
      temp_mat=nu_mat.rows(1,n1-1)%repmat((1-exp(ydev.subvec(1,n1-1))),1,nu_mat.n_cols);
      gr=gr+sum(temp_mat,0).t()+(0.5*(1-exp(ydev(0)))*nu_mat.row(0)).t()+
        (0.5*(1-exp(ydev(n1)))*nu_mat.row(n1)).t();
      
      
      // h
      big_mat=repmat((nu_mat.rows(1,n1-1)).t(),nbeta,1) % (matrix_big_mat(nu_mat.rows(1,n1-1),repmat(jj,nbeta,1)).t());
      
      coef_mat=repmat(exp(ydev.subvec(1,n1-1).t()),pow(nbeta,2),1);
      
      arma::mat temp_one=big_mat % coef_mat;
      arma::cube c_one(nbeta, nbeta, n1-1);
      
      int jj=0;
      
      for(int z=0;z<c_one.n_slices;z++)
      {
        for(int j=0;j<c_one.n_cols;j++)
        {
          for(int q=0;q<c_one.n_rows;q++)
          {
            c_one(q,j,z)=temp_one(jj);
            jj=jj+1;
            
          }
          
        }
      }
      
      
      arma::mat sum_mat_one=sum(c_one,2);
      he=he+sum_mat_one+
        (0.5*exp(ydev(0))*nu_mat.row(0)).t() * nu_mat.row(0)+
        (0.5*exp(ydev(n1))*nu_mat.row(n1)).t() * nu_mat.row(n1);
      
    }
    
  }
  
  
  
  
  Rcpp::List re = List::create(Named("value",0), _["gradient"]=0, _["hessian"]=0);
  
  
  re[0]=f;
  re[1]=gr;
  re[2]=he;
  
  return(re);
  
}
