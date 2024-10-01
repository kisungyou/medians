#include <RcppArmadillo.h>
#include "src_sphere_operations.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;




// =============================================================================
// auxiliary functions for spherical coordinates
// =============================================================================
arma::rowvec sphere_exp(arma::rowvec x, arma::rowvec d, double t){
  double nrm_td = arma::norm(t*d, 2);
  arma::rowvec output;
  if (nrm_td < 1e-12){
    output = x;
  } else{
    output = std::cos(nrm_td)*x + ((std::sin(nrm_td))/nrm_td)*t*d;
    output /= arma::norm(output, 2);
  }
  return(output);
}
arma::rowvec sphere_proj(arma::rowvec x, arma::rowvec u){
  return(u-x*(arma::dot(x,u)));
}
double sphere_dist(arma::rowvec x, arma::rowvec y){
  arma::rowvec vecxy = x-y;
  double dotxy = arma::dot(x, y);
  
  if (arma::norm(vecxy, 2) < arma::datum::eps*10.0){
    return(0.0);
  } else if (std::sqrt(dotxy*dotxy) >= (1.0 - 10.0*arma::datum::eps)){
    return(arma::datum::pi);
  } else {
    return(std::acos(arma::dot(x,y)));
  }
}
arma::rowvec sphere_log(arma::rowvec x, arma::rowvec y){
  arma::rowvec v = sphere_proj(x,y-x);
  double di = sphere_dist(x,y);
  if (di > 1e-6){
    double nv = arma::norm(v, 2);
    v = v*(di/nv);
  }
  return(v);
}
double sphere_metric(arma::rowvec x, arma::rowvec d1, arma::rowvec d2){
  return(arma::dot(d1,d2));
}
