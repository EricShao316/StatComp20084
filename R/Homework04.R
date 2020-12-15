#' @title importance functions
#' @name importfunctions
#' @description Homework04
#' @useDynLib StatComp20084
#' @examples
#' \dontrun{
#' set.seed(1207)
#' n = 1e4
#' # generate X from f1
#' x1 = rnorm(n)
#' # generate X from f2 (inverse transform method)
#' x2 = sqrt(-2*log(runif(n)))
#' # variance
#' var_f1 = var(x1^2)
#' var_f2 = var(x2)/(2*pi)
#' list(var_f1 = var_f1, var_f2 = var_f2)
#' g = function(x){
#' exp(-x - log(1+x^2)) * (x>0) * (x<1)
#' }
#' # f3 in Ex5.10 (inverse transform method)
#' x3 = -log(1 - runif(n) * (1 - exp(-1)))
#' f3g = g(x3) / (exp(-x3) / (1 - exp(-1)))
#' theta_hat_f3 = mean(f3g)
#' se_f3 = sd(f3g)
#' # Stratified importance sampling (inverse transform method)
#' k = 5
#' E = numeric(k)
#' var_sf = numeric(k)
#' gfi = numeric(n)
#' for(i in 1:k){
#'   y = runif(n, min = (i-1)/k, max = i/k)
#'   x = -log(exp((1-i)/k) - y * (exp((1-i)/k) - exp(-i/k)))
#'   gfi = g(x) / (exp(-x) / (exp((1-i)/k) - exp(-i/k)))
#'   E[i] = mean(gfi)
#'   var_sf[i] = var(gfi)
#' }
#' theta_hat_sf = sum(E)
#' se_sf = sqrt(sum(var_sf))
#' list(theta_hat_f3 = theta_hat_f3, theta_hat_sf = theta_hat_sf, 
#'      se_f3 = se_f3, se_sf = se_sf)
#' }
NULL