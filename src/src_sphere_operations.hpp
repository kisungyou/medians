#ifndef SRC_SPHERE_OPERATIONS_
#define SRC_SPHERE_OPERATIONS_

#include "RcppArmadillo.h"

arma::rowvec sphere_exp(arma::rowvec x, arma::rowvec d, double t);
arma::rowvec sphere_proj(arma::rowvec x, arma::rowvec u);
double sphere_dist(arma::rowvec x, arma::rowvec y);
arma::rowvec sphere_log(arma::rowvec x, arma::rowvec y);
double sphere_metric(arma::rowvec x, arma::rowvec d1, arma::rowvec d2);

#endif