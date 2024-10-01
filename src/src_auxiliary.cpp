#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// aux_row_normalize : row normalization to have unit norms


// =============================================================================
// aux_row_normalize
// =============================================================================
// [[Rcpp::export]]
arma::mat aux_row_normalize(arma::mat &X){
  int N = X.n_rows;
  int P = X.n_cols;
  
  arma::mat output(N,P,fill::zeros);
  arma::rowvec Xrow(P,fill::zeros);
  for (int n=0; n<N; n++){
    Xrow = X.row(n);
    output.row(n) = Xrow/arma::norm(Xrow, 2);
  }
  return(output);
}