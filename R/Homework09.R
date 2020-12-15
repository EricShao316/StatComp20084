#' @title eval_f0 Mle
#' @description Mle Calculation
#' @param x,x1,n.A,n.B,nOO,nAB Mle Constraction
#' @useDynLib StatComp20084
#' @export
eval_f0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}



#' @title eval_g0 Mle
#' @description Mle Calculation
#' @param x,x1,n.A,n.B,nOO,nAB Mle Constraction
#' @export
eval_g0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  return(sum(x)-0.999999)
}



#' @title sd.calculation
#' @description sd.calcu prepered to use function \code{Map}
#' @param x,x.bar using for sd.calculation
#' @export
sd.calcu = function(x,x.bar) {sum((x-x.bar)^2) / (length(x) - 1)
}