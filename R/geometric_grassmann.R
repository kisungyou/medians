#' Geometric Median on Grassmann Manifold
#' 
#' Grassmann manifold \eqn{Gr(p,k)} is a collection of \eqn{k}-dimensional subspaces in \eqn{\mathbb{R}^p}. 
#' For a collection of \eqn{k}-frames, it computes a weighted version of the geometric median 
#' under the chordal distance.
#' 
#' @param X X an \eqn{(p\times k\times n)} array that stacks \eqn{n} \eqn{k}-frames. It must be provided as a 3-dimensional array.
#' @param weights \code{NULL} for equal weight \code{rep(1/n,n)}; otherwise, it has to be a vector of length \eqn{n}.
#' @param ... extra parameters:\describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{Cauchy stopping criterion (default: 1e-8).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{solution}{an optimal solution as an \eqn{(p\times k)} orthonormal basis.}
#' \item{time}{elapsed wall-clock time for computation measured in seconds.}
#' }
#' 
#' @examples
#' \donttest{
#' ## load the iris data and extract the data
#' data("iris")
#' 
#' X = as.matrix(iris[,1:4])
#' lab = as.factor(iris[,5])
#' 
#' ## repeat 80% sampling and getting 2-dimensional PCA
#' proj3 <- array(0,c(4,2,50)) # repeat this 50 times
#' for (it in 1:50){
#'   # sampling
#'   Xsam = X[sample(1:nrow(X),round(nrow(X)*0.8)),]
#'   
#'   # compute eigenvector of empirical covariance
#'   Xcov = stats::cov(Xsam)
#'   Xeig = base::eigen(Xcov)
#'   
#'   # same the top-2 PCs
#'   proj3[,,it] = Xeig$vectors[,1:2]
#' }
#' 
#' ## compute geometric median
#' grmed = gmed_grassmann(proj3)
#' proj_gmed = grmed$solution
#' 
#' ## compute PCA for the original data
#' proj_base = base::eigen(stats::cov(X))$vectors[,1:2]
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(X%*%proj_base, col=lab, pch=19, xlab="", ylab="", main="base PCA")
#' plot(X%*%proj_gmed, col=lab, pch=19, xlab="", ylab="", main="geometric median")
#' par(opar)
#' }
#' 
#' @concept geometry
#' @export
gmed_grassmann <- function(X, weights=NULL, ...){
  ###############################################
  # Preprocessing
  if (!is.array(X)){
    stop("* gmed_grassmann: 'X' should be an array.")
  }
  if (length(dim(X))!=3){
    stop("* gmed_grassmann: 'X' should be a 3-dimensional array.")
  }
  n = dim(X)[3]

  # Parameters : explicit
  if ((length(weights)==0)&&(is.null(weights))){
    par_weights = rep(1/n, n)
  } else {
    if ((!is.vector(weights))||(length(weights)!=n)){
      stop(paste0("* gmed_grassmann : 'weights' should be a vector of length ",n))
    }
    if (any(weights <= 0)){
      stop(paste0("* gmed_grassmann : 'weights' should have nonnegative values."))
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
  t_init  = Sys.time()
  run_sol = cpp_geometric_grassmann(X, par_weights, par_maxiter, par_abstol)
  t_end   = Sys.time()
  duration = as.numeric(difftime(t_end, t_init, units="secs"))
  
  ###############################################
  # RETURN
  output = list(solution=as.matrix(run_sol),
                time=duration)
  return(output)
}
