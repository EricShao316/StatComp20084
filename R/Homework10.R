# Random Walk Metropolits sampler
#' @title Laplace density function R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{dlapR} and \code{rwlapR}) and Cpp functions (\code{gibbsC} and \code{vaccC}).
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm runif
#' @useDynLib StatComp20084
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark::microbenchmark(
#'   lapR = dlapR(10),
#'   lapC = dlapC(10)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' 
#' tm2 <- microbenchmark::microbenchmark(
#'   rwR = rwlapR(5000, 2),
#'   rwC = rwlapC(5000, 2)
#' )
#' print(summary(tm2)[,c(1,3,5,6)])
#' }
NULL

#' @title Laplace density function using R.
#' @description Laplace density function using R.
#' @param x the sample points (numeric)
#' @return the value of density function in x
#' @examples
#' \dontrun{
#' dla <- dlapR(10)
#' dla
#' }
#' @export
dlapR = function(x) .5 * exp(-abs(x))

#' @title A Random Walk sampler using R
#' @description A Random Walk sampler using R
#' @param N the length of chains
#' @param sigma the variance of normal distribution
#' @param x0 initail of x, default is 0
#' @return the generated chain
#' @examples
#' \dontrun{
#' rwlap <- rwlapC(5000, 2)
#' plot(1:5000, rwlap$x, ylab = "x", 
#'      main = expression(paste(sigma, '= 2')), type = "l")
#' }
#' @export
rwlapR = function(N, sigma, x0 = 0){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N){
    y = rnorm(1, x[i-1], sigma)
    if(u[i] <= ( dlapR(y) / dlapR(x[i-1]) )) x[i] = y
    else {x[i] = x[i-1] 
    k = k + 1}
  }
  return(list(x = x, k = k))
}