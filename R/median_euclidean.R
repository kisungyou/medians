#' Geometric Median in Euclidean Space
#' 
#' This function solves a weighted version of the L1-median problem, which is 
#' geometric median under the Euclidean setting. For a given data 
#' \eqn{x_1, x_2, \ldots, x_n \in \mathbb{R}^p}, it minimizes the following functional:
#' \deqn{F(\mu) = \sum_{i=1}^n w_i \|x_i - \mu\|_2,}
#' for some nonnegative weights \eqn{w_1, \ldots, w_n \in \mathbb{R}}.
#' 
#' @param X an \eqn{(n\times p)} matrix for \eqn{p}-dimensional signal. If vector is given, it is assumed that \eqn{p=1}.
#' @param weights \code{NULL} for equal weight \code{rep(1/n,n)}; otherwise, it has to be a vector of length \eqn{n}.
#' @param ... extra parameters:\describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{Cauchy stopping criterion (default: 1e-8).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{solution}{an optimal solution as a length-\eqn{p} vector.}
#' \item{time}{elapsed wall-clock time for computation measured in seconds.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## set seed for reproducibility
#' set.seed(496)
#' 
#' ## generate sin(x) data with noise for 100 replicates
#' X = seq(from=0,to=10,length.out=20)
#' Y = array(0,c(100,20))
#' for (i in 1:100){
#'    Y[i,] = sin(X) + stats::rnorm(20, sd=0.75)
#' }
#' 
#' ## compute L1-median and L2-mean
#' vecL2 = base::colMeans(Y)
#' vecL1 = median_euclidean(Y)$solution
#' 
#' ## compute deviation from the signal
#' signal = sin(X)
#' devL2  = sqrt(sum(vecL2-signal)^2)
#' devL1  = sqrt(sum(vecL1-signal)^2)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' matplot(t(Y[1:5,]), type="l", main="5 generated data", ylim=c(-2,2))
#' plot(t, vecL2, type="l", col="blue", main=paste0("L2:error=",round(devL2,4)), ylim=c(-2,2))
#' plot(t, vecL1, type="l", col="red",  main=paste0("L1:error=",round(devL1,4)), ylim=c(-2,2))
#' par(opar)
#' }
#' 
#' @export
median_euclidean <- function(X, weights=NULL, ...){
  ###############################################
  # Preprocessing
  if (is.vector(X)){
    X = matrix(X, ncol = 1)
  }
  n = nrow(X)
  d = ncol(X)
  
  # Parameters : explicit
  if ((length(weights)==0)&&(is.null(weights))){
    par_weights = rep(1/n, n)
  } else {
    if ((!is.vector(weights))||(length(weights)!=n)){
      stop(paste0("* median_euclidean : 'weights' should be a vector of length ",n))
    }
    if (any(weights) <= 0){
      stop(paste0("* median_euclidean : 'weights' should have nonnegative values."))
    }
    par_weights = weights/base::sum(weights);
  }
  
  # Parameters : implicit
  params  = list(...)
  pnames  = names(params)
  
  if ("maxiter" %in% pnames){
    par_maxiter = max(5, round(params$maxiter))
  } else {
    par_maxiter = 100
  }
  if ("abstol" %in% pnames){
    par_abstol = max(100*.Machine$double.eps, abs(params$abstol))
  } else {
    par_abstol = sqrt(.Machine$double.eps) # 1e-8
  }
  
  ###############################################
  # MAIN COMPUTATION & RETURN
  t_init   = Sys.time()
  run_sol  = cpp_geometric_euclidean(X, par_weights, par_maxiter, par_abstol)
  t_end    = Sys.time()
  duration = as.numeric(difftime(t_end, t_init, units="secs"))
    
  ###############################################
  # RETURN
  output = list(solution=as.vector(run_sol),
                time=duration)
  return(output)
}
