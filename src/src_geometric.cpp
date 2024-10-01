#include <RcppArmadillo.h>
#include "src_sphere_operations.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// cpp_geometric_euclidean
// cpp_geometric_grassmann  : just the chordal distance
// cpp_geometric_sphere_ext : spherical median with chordal distance
// cpp_geometric_sphere_int : spherical median with geodesic distance



// =============================================================================
// cpp_geometric_sphere_int
// =============================================================================
// [[Rcpp::export]]
arma::rowvec cpp_geometric_sphere_int(arma::mat &X, arma::vec weights, int maxiter, double abstol){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // prepare
  arma::rowvec Sold = arma::mean(X, 0);
  Sold /= arma::norm(Sold, 2);
  arma::rowvec Snew(P,fill::zeros);
  double Sinc = 0.0;
  
  arma::rowvec Stmp(P,fill::zeros);
  arma::mat Slogs(N,P,fill::zeros);
  arma::vec Sdist(N,fill::zeros);
  
  arma::rowvec tmp1(P,fill::zeros);
  double tmp2 = 0.0;
  
  double epsnum = arma::datum::eps*100.0; // machine epsilon; approx. 2*1e-14
  int M = 0;
  
  // iteration
  arma::uvec nonsingular;
  // main iteration
  for (int it=0; it<maxiter; it++){
    // 1. compute log-pulled vectors and norm
    Slogs.fill(0.0);
    Sdist.fill(0.0);
    for (int n=0; n<N; n++){
      Stmp = sphere_log(Sold, X.row(n));
      Slogs.row(n) = Stmp;
      Sdist(n) = std::sqrt(sphere_metric(Sold, Stmp, Stmp));
    }
    
    // 2. find the one with singular-distance
    nonsingular.reset();
    nonsingular = arma::find(Sdist > epsnum);
    M = nonsingular.n_elem;
    if (M < 1){
      break;
    }
    
    // 3. update numerator and denominator
    tmp1.fill(0.0);
    tmp2 = 0.0;
    for (int j=0; j<M; j++){
      tmp1 += weights(nonsingular(j))*Slogs.row(nonsingular(j))/Sdist(nonsingular(j));
      tmp2 += weights(nonsingular(j))/Sdist(nonsingular(j));
    }
    
    // 4. update to the new one
    Stmp = tmp1/tmp2;
    Snew = sphere_exp(Sold, Stmp, 1.0);
    Sinc = arma::norm(Snew-Sold, "fro");
    Sold = Snew;
    
    // 5. update information
    if (Sinc < abstol){
      break;
    }
  }
  
  // return
  return(Sold);
}



// =============================================================================
// cpp_geometric_sphere_ext 
// =============================================================================
// [[Rcpp::export]]
arma::rowvec cpp_geometric_sphere_ext(arma::mat &X, arma::vec weights, int maxiter, double abstol){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // prepare
  arma::rowvec x_old = arma::mean(X, 0);
  x_old /= arma::norm(x_old, 2);
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
    x_new /= arma::norm(x_new, 2);
    
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


// =============================================================================
// cpp_geometric_euclidean 
// =============================================================================
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




// =============================================================================
// cpp_geometric_grassmann : just the chordal distance 
// =============================================================================
// grass_eigenK         : compute eigendecomposition and select first K vectors
// grass_genQ           : generate a (p \times k) orthogonal matrix
// cpp_median_grassmann : given an (p \times k \times n) array, use extrinsic.

arma::mat grass_eigenK(arma::mat P, int K){
  arma::vec eigval;
  arma::mat eigvec;
  
  /// eigenvalue decomposition in an ascending order
  arma::eig_sym(eigval, eigvec, P);
  
  // hence, we must choose the last columns
  arma::mat output = eigvec.tail_cols(K);
  return(output);
}
arma::mat grass_genQ(int p, int k){
  // create a random matrix
  arma::mat X(p,k,arma::fill::randn);
  arma::mat Q; 
  arma::mat R;
  
  arma::qr(Q,R,X);
  return(Q);
}

// [[Rcpp::export]]
arma::mat cpp_geometric_grassmann(arma::cube &X, arma::vec weights, int maxiter, double abstol){
  // parameters
  int P = X.n_rows;
  int K = X.n_cols;
  int N = X.n_slices;
  
  // prepare
  arma::mat Q_old = grass_genQ(P,K);
  arma::mat P_old = Q_old*Q_old.t();
  arma::mat Q_new(P,K,fill::zeros);
  arma::mat P_new(P,P,fill::zeros);
  double P_inc = 0.0;
  
  arma::vec dists(N,fill::zeros);
  arma::mat tmp_P(P,P,fill::zeros);
  
  double tmp_c = 0.0;
  double epsnum = arma::datum::eps*100.0; // machine epsilon; approx. 2*1e-14
  
  // main iteration
  for (int it=0; it<maxiter; it++){
    // step 1. compute distances; if too close to anchors, stop.
    for (int n=0; n<N; n++){
      dists(n) = arma::norm(P_old - X.slice(n)*arma::trans(X.slice(n)), "fro");
    }
    if (dists.min() < epsnum){
      break;
    }
    
    // step 2. compute numerator and denominator
    tmp_P.fill(0.0);
    tmp_c = 0.0;
    for (int n=0; n<N; n++){
      tmp_c += weights(n)/dists(n);
      tmp_P += (weights(n)/dists(n))*(X.slice(n)*arma::trans(X.slice(n)));
    }
    P_new = (tmp_P/tmp_c);
    Q_new = grass_eigenK(P_new, K);
    
    // step 3. updating information
    P_inc = arma::norm(P_new - P_old, "fro");
    P_old = P_new;
    Q_old = Q_new;
    if (P_inc < abstol){
      break;
    }
  }
  
  // return
  return(Q_old);
}