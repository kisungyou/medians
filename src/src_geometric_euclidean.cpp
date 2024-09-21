#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::rowvec cpp_geometric_euclidean(arma::mat &X, arma::vec weights, int maxiter, double abstol){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // prepare
  arma::rowvec x_old = arma::mean(X, 0);
  arma::rowvec x_new(P,fill::zeros);
  double x_inc = 0.0;
  
  arma::vec dists(N,fill::zeros);
  arma::rowvec tmp_x(P,fill::zeros);
  double tmp_c = 0.0;
  double epsnum = arma::datum::eps*100.0; // machine epsilon; approx. 2*1e-14
  
  // main iteration
  for (int it=0; it<maxiter; it++){
    // step 1. compute distances; if too close to one of the points, stop.
    for (int n=0; n<N; n++){
      dists(n) = arma::norm(X.row(n) - x_old, 2);
    }
    if (dists.min() < epsnum){
      break;
    }
    
    // step 2. compute numerator and denominator
    tmp_x.fill(0.0);
    tmp_c = 0.0;
    for (int n=0; n<N; n++){
      tmp_c += weights(n)/dists(n);
      tmp_x += (weights(n)/dists(n))*X.row(n);
    }
    x_new = tmp_x/tmp_c;
    
    // step 3. updating information
    x_inc = arma::norm(x_new-x_old, 2);
    x_old = x_new;
    if (x_inc < abstol){
      break;
    }
  }
  
  // return
  return(x_old);
}
