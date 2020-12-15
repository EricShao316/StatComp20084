#' @title MC estimate
#' @name MC
#' @description Homework03
#' @useDynLib StatComp20084
#' @examples
#' \dontrun{
#' n1 = 1e5
#' set.seed(2320)
#' u1 = runif(n1, min = 0, max = pi/3)
#' theta1 = sum(pi/3*sin(u1))/n1
#' cp1 = c(theta1, 0.5)
#' cp1
#' n2 = 1e5
#' m = 1000
#' theta_mc = numeric(0)
#' theta_av = numeric(0)
#' per = numeric(0)
#' for (i in 1:m){
#'   u2 = runif(n2/2)
#'   u3 = 1 - u2
#'   theta_mc[i] = sum(exp(u2))/n2
#'   theta_av[i] = (sum(exp(u2)) + sum(exp(u3)))/n2
#' }
#' per = (var(theta_mc) - var(theta_av))/var(theta_mc)
#' exact_per = -(exp(1)-(exp(1)-1)^2)/0.2420357
#' list(var_mc = var(theta_mc), var_av = var(theta_av), per = per, exact_per = exact_per )
#' }
NULL