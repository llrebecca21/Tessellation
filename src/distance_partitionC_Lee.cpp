#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec distance_partitionC_Lee(arma::mat X, arma::vec centers, arma::vec weights)
{
  // define variables
  // define number rows of X
  int n = X.n_rows;
  // define number elements in centers
  int M = centers.n_elem;
  // define number of columns in X
  int p = X.n_cols;
  // initialize indices
  int ind;
  
  // Create empty vectors
  arma::vec prt(n, arma::fill::zeros);
  arma::vec dists(M, arma::fill::zeros);
  arma::uvec sort(M, arma::fill::zeros);
  // initialize the zero vector p with zeros
  arma::vec distance(p, arma::fill::zeros);

  
  if(true){
    arma::mat difference(n,p);
    arma::rowvec weights2 = (weights % weights).t();
    arma::vec current_dist(n);
    arma::vec smallest_dist(n , arma::fill::value(1e12) );
    for(int jj = 0 ; jj < M ; jj++){
      //calculate difference of X row with center_jj
      difference = (X.each_row() - X.row(centers(jj) - 1));
      // square the difference matrix
      difference %= difference;
      // calculate weights^2 * difference
      difference.each_row() %= weights2;
      // sum across the rows of the difference matrix
      current_dist = sum(difference, 1);
      // create uvec to identify new closest centers
      arma::uvec newclose = find(current_dist < smallest_dist);
      // make changes to prt at these indices
      prt(newclose) = arma::vec(newclose.n_elem , arma::fill::value(jj + 1));
      // make changes to smallest_dist at these indices
      smallest_dist(newclose) = current_dist(newclose);
    }
    return(prt);
  }
  
  if(false){
    // first for loop: each row in X
    for(int ii = 0 ; ii < n ; ii++)
    {
      // second for loop: for each center 
      for(int jj = 0 ; jj < M ; jj++)
      {
        // third for loop: for each column in X
        for(int zz = 0 ; zz < p ; zz++)
        {
          // Calculate Euclidean distance: element-by-element
          distance(zz)=pow((X(ii,zz)-X(centers(jj)-1,zz)),2);
        }
        
        // Create full weighted Euclidean distance between the two vectors (row of X_ii and center_jj)
        dists(jj)=sum(pow(weights,2)%distance);
        
      }
      
      // sort the dists
      sort=sort_index(dists);
      // choose first sorted index
      ind=sort(0);
      // reindex ind for R
      prt(ii)=ind+1;
      
    }
    
    //return(distance);
    //return the indices
    return(prt);
  }
}







