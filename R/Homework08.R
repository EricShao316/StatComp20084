#' @title Density of Laplace distribution
#' @name dlap
#' @param x vector (numeric)
#' @return the density value at x
#' @export
dlap = function(x) .5 * exp(-abs(x))



#' @title Random Walk Metropolits sampler
#' @description Random Walk Metropolits sampler for laplace distribution
#' @param N length of chain
#' @param sigma variance of normal distribution
#' @param x0 initail of x, default value is 0
#' @return the generated chain
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' n = 2000
#' rl1 = rw.lap(N = n, sigma = .05, x0 = 0)
#' rl2 = rw.lap(N = n, sigma = .5, x0 = 0)
#' rl3 = rw.lap(N = n, sigma = 2, x0 = 0)
#' rl4 = rw.lap(N = n, sigma = 16, x0 = 0)
#' knitr::kable(data.frame(k1 = rl1$k, k2 = rl2$k, k3 = rl3$k, k4 = rl4$k))
#' par(mfrow=c(2,2))
#' plot(1:n, rl1$x, ylab = "x", main = expression(paste(sigma, '= 0.05')), type = "l")
#' plot(1:n, rl2$x, ylab = "x", main = expression(paste(sigma, '= 0.5')), type = "l")
#' plot(1:n, rl3$x, ylab = "x", main = expression(paste(sigma, '= 2')), type = "l")
#' plot(1:n, rl4$x, ylab = "x", main = expression(paste(sigma, '= 16')), type = "l")
#' }
#' @export
rw.lap = function(N, sigma, x0 = 0){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N){
    y = rnorm(1, x[i-1], sigma)
    if(u[i] <= ( dlap(y) / dlap(x[i-1]) )) x[i] = y
    else {x[i] = x[i-1] 
    k = k + 1}
  }
  return(list(x = x, k = k))
}



#' @title Function to compute the diagnnostic statistic R
#' @description compute the diagnnostic statistic R according to Gelman-Rubin method.
#' @param psi the woking matrix used to calculate R.hat
#' @importFrom stats var
#' @return R.hat
#' @export
GR = function(psi){
  psi = as.matrix(psi)
  psi.si = numeric(k)
  n = ncol(psi)
  k = nrow(psi)
  psi.means = rowMeans(psi)
  B = n * var(psi.means)
  psi.w = apply(psi, 1, "var")
  W = mean(psi.w)
  Var.psi = (n-1) / n * W + B / n
  R.hat = Var.psi / W
  return(R.hat)
}



#' @title M-H sampler 
#' @description M-H sampler 
#' @importFrom stats rnorm runif dnorm
#' @param N length of chain
#' @param sigma variance of normal distribution
#' @param x0 initail of x, default value is 0
#' @examples 
#' \dontrun{
#' # Choose initial values
#' k = 4
#' X0 = c(-10, -5, 5, 10)
#' n = 15000
#' sigma = 2
#' 
#' # Generate the chains
#' X = matrix(NA, k, n)
#' for(i in 1:k) X[i,] = rw.lap(N = n, sigma = sigma, x0 = X0[i])$x
#' 
#' # Compute the diagnostic statistics: mean of ith chain up to time j
#' psi = t(apply(X, 1, "cumsum"))
#' for(i in 1:k) psi[i,] = psi[i,] / (1 : ncol(psi))
#' R.hat = GR(psi)
#' print(R.hat)
#' }
#' \dontrun{
#' # Plot of Chains psi
#' par(mfrow = c(2,2))
#' for(i in 1:k){
#'   plot(psi[i,1000:n],type = "l",ylab = expression(paste(psi)) )
#' }
#' #Plot of R.hat statistics
#' rhat = numeric(n)
#' for(j in 501:n) rhat[j] = GR(psi[,1:j])
#' par(mfrow = c(1,1))
#' plot(rhat[501:n], type = "l", xlab = " ", ylab = "R")
#' abline(h = 1.2, lty = 4)
#' }
#' @export
mh.lap = function(N, sigma, x0 = 0){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N){
    y = rnorm(1, x[i-1], sigma)
    if(u[i] <= ( (dlap(y) * dnorm(x[i-1], y, sigma)) / (dlap(x[i-1]) * dnorm(y, x[i-1], sigma)) )) x[i] = y
    else {x[i] = x[i-1] 
    k = k + 1}
  }
  return(list(x = x, k = k))
}



#' @title Root finder
#' @description Root finding function
#' @importFrom stats uniroot
#' @param k parameter vector k
#' @return unitroot at parameter k
#' @examples 
#' \dontrun{
#' set.seed(2156)
#' k = c(4:25, 100, 500, 1000)
#' kS = data.frame(k = k , S.root = S.root(k = k))
#' par(mfrow=c(2,1))
#' plot(kS$k[1:22], kS$S.root[1:22], type = "l", ylab = "S.root", xlab = "k", lwd = 2, col = "darkred")
#' plot(kS$k, kS$S.root, type = "l", ylab = "S.root", xlab = "k", lwd = 2, col = "darkred")
#' knitr::kable(kS[-(3:18),])
#' }
#' @export
S.root = function(k){
  n = length(k)
  sroot = numeric(n)
  for(i in 1:n){
    f = function(x) S1(x,k[i]) - S2(x,k[i])
    sroot[i] = uniroot(f,c(1e-5,sqrt(k[i])-1e-5))$root
  }
  sroot
}



#' @title function of the first curve
#' @importFrom stats pt
#' @param a,k parameter a,k
#' @export
S1 = function(a, k) 1 - pt(sqrt(a^2*(k-1)/(k-a^2)), df = k-1)



#' @title function of the second curve
#' @importFrom stats pt
#' @param a,k parameter a,k
#' @export
S2 = function(a, k) 1 - pt(sqrt(a^2*k/(k+1-a^2)), df = k)



