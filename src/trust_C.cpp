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
Rcpp::List trust_C(arma::vec parinit, double rinit, double rmax, arma::vec parscale, int iterlim, double fterm, double mterm, bool minimize, arma::vec index, arma::vec nseg_time_temp, Rcpp::List ytemp, double tau_temp,int nbeta,int nbasis, double sigmasqalpha)
{
  
  double r=rinit;
  bool rescale=true;
  arma::vec theta=parinit;
  
  //Rcpp::List out=whittle_derivs2_C(parinit,index, nseg_time_temp, ytemp, tau_temp, nbeta, nbasis, sigmasqalpha);
  Rcpp::List out=whittle_derivs2_C(theta,index, nseg_time_temp, ytemp, tau_temp, nbeta, nbasis, sigmasqalpha);
  
  bool accept=true;
  
  
  int iiter;
  double f=0.0;
  arma::vec g;
  arma::mat B;
  
  
  NumericMatrix B_mat;
  NumericVector g_vec;
  arma::mat Eivector;
  arma::vec Eivalue;
  arma::mat gq;
  arma::vec ptry;
  arma::vec wtry;
  arma::vec utry;
  arma::vec vtry;
  arma::vec theta_try;
  
  bool is_hard;
  bool is_easy;
  bool is_newton;
  bool is_terminate;
  double rho=0;
  
  double C1;
  double C2;
  double C3;
  double uout;
  double ftry;
  
  for(iiter=0;iiter<iterlim;iiter++)
  {
    
    
    if(accept)
    {
      // output B
      NumericVector B_vec=out["hessian"];
      B_vec.attr("dim") = Dimension(nbeta, nbeta);
      B_mat=as<NumericMatrix>(B_vec);
      B=as<arma::mat>(B_mat);
      
      // output g
      g_vec=out["gradient"];
      g=as<arma::vec>(g_vec);
      
      // output f
      f=out["value"];
      
      if(rescale)
      {
        B=B/rcpparma_outerproduct(parscale);
        g=g/parscale;
      }
      
      if (!minimize) {
        
        B=-B;
        g=-g;
        f=-f;
        
      }
      
      
      Eivector=eigenVector(B);
      Eivalue=eigenValue(B);
      
      gq=Eivector.t()*g;
      
      
    }
    
    
    // solve trust region subproblem
    // try for Newton
    is_newton=false;
    if (all(Eivalue>0))
    {
      
      ptry=-(Eivector*(gq/Eivalue));
      
      
      if(norm(ptry)<=r)
      {
        is_newton=true;
      }
      
      
    }
    
    
    // non-Newton
    if(!is_newton)
    {
      
      
      double lambda_min=min(Eivalue);
      
      arma::vec trust_beta=Eivalue-lambda_min;
      
      
      arma::uvec inotmin= find(trust_beta!=0);
      arma::uvec imin= find(trust_beta==0);
      
      arma::vec ratio=gq/trust_beta;
      
      C1=sum(pow(ratio(inotmin),2));
      C2=sum(pow(gq(imin),2));
      arma::vec temp=sum(pow(gq,2));
      C3=temp(0);
      
      
      
      if( (C2>0) || (C1>pow(r,2)) )
      {
        
        
        is_easy= true;
        is_hard= (C2==0);
        
        // easy cases
        double beta_dn=std::sqrt(C2)/r;
        double beta_up=std::sqrt(C3)/r;
        
        
        
        if(fredC(beta_up,gq,trust_beta,C1,C2,r)<=0)
        {
          uout=beta_up;
          
        }
        else if(fredC(beta_dn,gq,trust_beta,C1,C2,r)>=0)
        {
          uout=beta_dn;
          
        }else{
          
          
          Rcpp::List root_result=unirootC(beta_dn,beta_up,gq,trust_beta,C1,C2,r);
          uout=root_result["root"];
          
          
        }
        
        wtry=gq/(trust_beta+uout);
        ptry=-Eivector*wtry;
        
        
        
      }else{
        
        is_hard=true;
        is_easy=false;
        
        // hard-hard case
        wtry=gq/trust_beta;
        wtry.elem(imin).fill(0);
        ptry= -Eivector*wtry;
        utry=sqrt(pow(r,2)-sum(pow(ptry,2)));
        
        if(utry(0)>0){
          
          vtry=Eivector.col(imin(0));
          ptry=ptry+utry*vtry;
          
        }
        
      }
      
    }
    
    
    // predicted versus actual change
    double preddiff = sum(ptry % (g + (B * ptry) / 2));
    
    if (rescale) {
      theta_try = theta + ptry / parscale;
    } else {
      theta_try = theta + ptry;
    }
    
    
    out=whittle_derivs2_C(theta_try,index, nseg_time_temp, ytemp, tau_temp, nbeta, nbasis, sigmasqalpha);
    ftry=out["value"];
    
    
    
    if (! minimize)
    {
      ftry = (- ftry);
    }
    
    
    rho = (ftry - f) / preddiff;
    
    
    // termination test
    if (ftry < arma::math::inf()) {
      is_terminate = (abs(ftry - f) < fterm) || (abs(preddiff) < mterm);
    } else {
      is_terminate = false;
      rho = (-1* arma::math::inf());
    }
    
    // adjustments
    if (is_terminate) {
      if (ftry < f) {
        accept = true;
        theta = theta_try;
      }
    } else {
      
      
      if (rho <  (1.0 / 4)) {
        accept = false;
        r = r / 4.0;
      } else {
        accept = true;
        theta = theta_try;
        if (rho > (3.0 / 4) && (! is_newton))
        {
          r = std::min(2*r, rmax);
        }
        
      }
      
    }
    
    if (is_terminate)
    {
      break;
    }
    
    
  }
  
  
  Rcpp::List out_re=whittle_derivs2_C(theta, index, nseg_time_temp, ytemp, tau_temp, nbeta, nbasis, sigmasqalpha);
  
  Rcpp::List re =List::create(Named("value",0), _["gradient"]=0, _["hessian"]=0,_["argument"]=0,_["converged"]=0, _["iterations"]=0);
  
  re[0]=out_re[0];
  re[1]=out_re[1];
  re[2]=out_re[2];
  re[3]=theta;
  re[4]=is_terminate;
  re[5]=iiter+1;
  
  
  return(re);
  
  
  
  
}
