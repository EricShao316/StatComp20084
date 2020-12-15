#' @title Different density plots of different sample size.
#' @name latticePlot
#' @description Different density plots of different sample size.
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @useDynLib StatComp20084
#' @examples
#' \dontrun{
#' n = seq(5, 45, 5)
#' x = rnorm(sum(n))
#' y = factor(rep(n, n), labels = paste("n =", n))
#' lattice::densityplot(~ x | y,
#' panel = function(x, ...) {
#' lattice::panel.densityplot(x, col = "DarkOliveGreen", ...)
#' lattice::panel.mathdensity(dmath = dnorm,
#'                    args = list(mean = mean(x), sd=sd(x)),
#'                    col = "darkblue")
#' })
#' }
NULL