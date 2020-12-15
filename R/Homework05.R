#' @title skewness calculation
#' @name skewness.calcul
#' @description skewness calculation
#' @useDynLib StatComp20084
#' @param x Sample data (numeric vector)
#' @return skewness
sk = function(x){
  xbar = mean(x)
  m3 = mean((x-xbar)^3)
  m2 = mean((x-xbar)^2)
  return(m3 / m2^(1.5))
}



#' @title sktest power calculation
#' @description sktest power calculation
#' @importFrom stats rnorm rbeta rt
#' @param n Sample size (integer)
#' @param m Times of replicate (integer)
#' @param cv Crit.values for each n
#' @param type The type of distribution.
#' @return mpowerx, power of sktest over m replicate.
#' @examples 
#' \dontrun{
#' # empirical power
#' m = 1000
#' alpha = 0.05
#' n = c(10, 20, 30, 50, 100, 500) # sample sizes
#' cv = qnorm(1-alpha, 0, sqrt(6/n)) #crit.values for each n
#' sknorm = sktest(n = n, m = m,cv = cv, type = "norm")
#' skbeta = sktest(n = n, m = m,cv = cv, type = "beta")
#' skt = sktest(n = n, m = m,cv = cv, type = "t")
#' list(sknorm = sknorm, skbeta = skbeta, skt = skt )
#' }
#' @export
sktest = function(n, m, cv, type = "norm"){
  mpowerx = numeric(length(n))
  if(type == "norm"){
    for(i in 1:length(n)){
      powerx = numeric(m)
      for(j in 1:m){
        x = rnorm(n[i], 2, 3)
        skx = sk(x)
        powerx[j] = abs(skx) > cv[i] 
      }
      mpowerx[i] = mean(powerx)
    }
    return(mpowerx)
  }
  if(type == "beta"){
    for(i in 1:length(n)){
      powerx = numeric(m)
      for(j in 1:m){
        x = rbeta(n[i], 2, 2)
        skx = sk(x)
        powerx[j] = abs(skx) > cv[i] 
      }
      mpowerx[i] = mean(powerx)
    }
    return(mpowerx)
  }
  if(type == "t"){
    for(i in 1:length(n)){
      powerx = numeric(m)
      for(j in 1:m){
        x = rt(n[i], df = 10)
        skx = sk(x)
        powerx[j] = abs(skx) > cv[i] 
      }
      mpowerx[i] = mean(powerx)
    }
    return(mpowerx)
  }
}



#' @title Count-5-test
#' @description The equal variance test by count5test.
#' @param x,y Two group of data in the same size.
#' @return Test statistic for count-5-test.
#' @examples 
#' \dontrun{
#' sigma1 = 1
#' sigma2 = 1.5
#' m = 1000
#' n = c(10, 100, 1000)
#' power_ct = numeric(length(n))
#' power_f = numeric(length(n))
#' for(i in 1:length(n)){
#'   power_ct[i] = mean(replicate(m, expr = {
#'     x = rnorm(n[i], 0, sigma1)
#'     y = rnorm(n[i], 0, sigma2)
#'     count5test(x, y)
#'   }))
#'   power_f[i] = mean(replicate(m, expr = {
#'     x = rnorm(n[i], 0, sigma1)
#'     y = rnorm(n[i], 0, sigma2)
#'     var.test(x, y)$p.value < 0.05
#'   }))
#' }
#' list(power_ct = power_ct, power_f = power_f)
#' }
#' @export
count5test = function(x, y){
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}



#' @title Mardia’s multivariate skewness test
#' @description  Mardia’s multivariate skewness test and compute type 1 error.
#' @param x Random Matrix
#' @importFrom stats cov
#' @return b1,d
#' @examples 
#' \dontrun{
#' m = 100
#' alpha = 0.05
#' n = c(10, 100, 1000)
#' d = 2 
#' mum = c(0, 0)
#' sigmam = matrix(c(1, 0.3, 0.3, 2), nrow = 2)
#' # type 1 error
#' cvm = qchisq(1 - alpha, d * (d+1) * (d+2) / 6)
#' p.reject = numeric(length(n))
#' for(i in 1:length(n)){
#'   skmtest = numeric(m)
#'   for(j in 1:m){
#'     x = MASS::mvrnorm(n[i], mum, sigmam)
#'     skmtest[j] = as.integer((skm(x) * n[i] / 6 )> cvm)
#'   }
#'   p.reject[i] = mean(skmtest)
#' }
#' p.reject
#' }
#' @export
skm = function(x){
  n = nrow(x)
  xbar = as.vector(apply(x, 2, mean))
  xbar = matrix(xbar, nrow = n, ncol = ncol(x))
  sigma = cov(x)
  insigma = solve(sigma)
  s = ((x - xbar) %*% insigma %*% t(x - xbar)) ^ 3
  return(mean(s))
}



#' @title Calculate the power of Mardia’s multivariate skewness test.
#' @name multivariate.skewness.test
#' @description the power of the skewness test
#' @importFrom stats qchisq
#' @examples 
#' \dontrun{
#' # power of skewness test
#' alpha <- 0.1 
#' n <- 30 
#' m <- 2500
#' d <- 2
#' epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05)) 
#' N <- length(epsilon) 
#' pwr <- numeric(N) 
#' cv <- qchisq(1-alpha,d * (d+1) * (d+2) / 6) 
#' for (j in 1:N) {
#'   #for each epsilon
#'   e <- epsilon[j]
#'   sktests <- numeric(m)
#'   for (i in 1:m) {
#'     #for each replicate
#'     x1 <- MASS::mvrnorm(n, c(0,0), matrix(c(1,0.3,0.3,1.5), nrow = 2)) 
#'     x2 <- MASS::mvrnorm(n, c(0,0), matrix(c(100,0.3,0.3,150), nrow = 2)) 
#'     Iprob = sample(c(1,2), n, replace = TRUE, prob = c(1 - e, e ))
#'     x = as.integer(Iprob == 1) * x1 + as.integer(Iprob == 2) * x2
#'     sktests[i] <- as.integer(n*skm(x)/6 >= cv) }
#'   pwr[j] <- mean(sktests) }
#' #plot power vs epsilon
#' plot(epsilon, pwr, type = "b",
#'      xlab = bquote(epsilon))
#' abline(h = .1, lty = 3)
#' se <- sqrt(pwr * (1-pwr) / m) 
#' #add standard errors
#' lines(epsilon, pwr+se, lty = 3)
#' lines(epsilon, pwr-se, lty = 3)
#' }
NULL


