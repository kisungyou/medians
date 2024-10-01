#' Geometric Median on Unit Hypersphere
#' 
#' Unit hypersphere \eqn{\mathbb{S}^{p-1}} is a collection of unit-norm vectors in 
#' \eqn{\mathbb{R}^p}, whose dimension is \eqn{p-1}. For a collection of \eqn{p}-dimensional 
#' unit norm vectors, it computes a weighted version of the geometric median 
#' under geodesic or chordal distance.
#' 
#' 
#' @param X an \eqn{(n\times p)} matrix whose rows are unit-norm vectors in \eqn{\mathbb{R}^p}.
#' @param weights \code{NULL} for equal weight \code{rep(1/n,n)}; otherwise, it has to be a vector of length \eqn{n}.
#' @param ... extra parameters:\describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{Cauchy stopping criterion (default: 1e-8).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{solution}{an optimal solution as a length-\eqn{p} unit-norm vector.}
#' \item{time}{elapsed wall-clock time for computation measured in seconds.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## set seed for reproducibility
#' set.seed(496)
#' 
#' ## generate sin(x) data with noise for 100 replicates
#' ## this time, each signal is normalized to have norm 1
#' X = seq(from=0,to=10,length.out=20)
#' Y = array(0,c(100,20))
#' for (i in 1:100){
#'    Y[i,] = sin(X) + stats::rnorm(20, sd=0.25)
#'    Y[i,] = Y[i,]/sqrt(sum(Y[i,]^2))
#' }
#' 
#' ## compute mean as a baseline for comparison
#' vecL2 = base::colMeans(Y)
#' 
#' ## compute medians under two geometries
#' vecLI = gmed_sphere(Y, geometry="intrinsic")$solution 
#' vecLE = gmed_sphere(Y, geometry="extrinsic")$solution
#' 
#' ## compute deviation from the normalized signal
#' signal = sin(X)
#' signal = signal/sqrt(sum(signal^2))
#' 
#' devL2  = sqrt(sum(vecL2-signal)^2)
#' devLE  = sqrt(sum(vecLE-signal)^2)
#' devLI  = sqrt(sum(vecLI-signal)^2)
#' 
#' ## visualize
#' yvis <- c(-0.5, 0.5) # custom range for y-axis plotting
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' matplot(t(Y[1:5,]), type="l", main="5 generated data", ylim=yvis)
#' plot(X, vecL2, type="l", col="green", main=paste0("mean:error=",round(devL2,4)), ylim=yvis)
#' lines(X, signal, lwd=1, lty=3)
#' plot(X, vecLE, type="l", col="red",  main=paste0("extrinsic median:error=",round(devLE,4)), ylim=yvis)
#' lines(X, signal, lwd=1, lty=3)
#' plot(X, vecLI, type="l", col="blue",  main=paste0("intrinsic median:error=",round(devLI,4)), ylim=yvis)
#' lines(X, signal, lwd=1, lty=3)
#' par(opar)
#' }
#' 
#' @concept geometric
#' @export
gmed_sphere <- function(X, weights=NULL, ...){
  ###############################################
  # Preprocessing
  if (!is.matrix(X)){
    stop("* gmed_sphere : an input 'X' must be a matrix.")
  }
  n = nrow(X)
  d = ncol(X)
  
  if (any(abs(rowSums(X^2)-rep(1,n)) > 10*.Machine$double.eps)){
    warning("* gmed_sphere : all rows must have unit norms. Normalization is applied.")
    X = aux_row_normalize(X)
  }
  
  # Parameters : explicit
  if ((length(weights)<1)&&(is.null(weights))){
    par_weights = rep(1/n, n)
  } else {
    if ((!is.vector(weights))||(length(weights)!=n)){
      stop(paste0("* gmed_sphere : 'weights' should be a vector of length ",n))
    }
    if (min(weights)<=0){
      stop(paste0("* gmed_sphere : 'weights' should have nonnegative values."))
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
  if ("geometry" %in% pnames){
    par_geometry = match.arg(tolower(params$geometry), c("intrinsic","extrinsic"))
  } else {
    par_geometry = "intrinsic"
  }
  
  ##############################################
  # COMPUTATION
  # case-branching run with time measurement
  t_init   = Sys.time()
  if (par_geometry=="intrinsic"){
    result = cpp_geometric_sphere_ext(X, par_weights, par_maxiter, par_abstol)
  } else {
    result = cpp_geometric_sphere_int(X, par_weights, par_maxiter, par_abstol)
  }
  t_end    = Sys.time()
  duration = as.numeric(difftime(t_end, t_init, units="secs"))


  ###############################################
  # RETURN
  output = list(solution=as.vector(result),
                time=duration)
  return(output)
}