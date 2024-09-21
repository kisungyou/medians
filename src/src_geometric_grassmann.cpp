#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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