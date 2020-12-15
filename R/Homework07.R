#' @title Comparsion between count5test and permutation version.
#' @name count5test.permutation
#' @description A Comparsion between permutation test \code{permu.count5} and ordinary count5test \code{count5test} for equal variance based on the maximum number of extreme points that applies when sample sizes are not necessarily equal.
#' @useDynLib StatComp20084
#' @importFrom stats rnorm
#' @examples 
#' \dontrun{
#' # MC-permutation-count 5 test
#' set.seed(1630)
#' n1 = c(30, 40, 45, 50, 105)
#' n2 = c(20, 40, 40, 40, 90)
#' mu1 = mu2 = 2
#' sigma1 = sigma2 = 1
#' m = 1000
#' permu.compare(group1 = n1, group2 = n2)
#' }
NULL


#' @title count 5 test for equal variance
#' @description ordinary count 5 test.
#' @param x,y two sample sets
#' @export
count5test = function(x, y){
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}



#' @title permutation-count 5 test for x,y
#' @description permutation version of count 5 test.
#' @param R times of replicate
#' @param n1,n2 two sample sizes
#' @param mu1,mu2 means of sample sets
#' @param sigma1,sigma2 true variance of two sample 
#' @export
permu.count5 = function(R = 999, n1 = 20, n2 = 30,
                        mu1 = 0, mu2 = 0,
                        sigma1 = 1, sigma2 = 1){
  x = rnorm(n1, mu1, sigma1)
  y = rnorm(n2, mu2, sigma2)
  z = c(x, y)
  count.hat = count5test(x, y)
  count.z = numeric(R)
  for(i in 1:R){
    k = sample(1 : (n1 + n2), n1, replace = FALSE)
    count.z[i] = count5test(z[k], z[-k])
  }
  alphahat = mean(count.z)
  return(alphahat)
}



#' @title comparsion between two methods
#' @description comparsion between two methods
#' @param group1,group2 two groups of samples
#' @return a data.frame for comparision
#' @export
permu.compare = function(group1, group2){
  mu1 = mu2 = 2
  sigma1 = sigma2 = 1
  m = 1000
  n = length(group1)
  alphahat.permu = numeric(n) 
  alphahat.no_permu = numeric(n)
  for(i in 1:n){
    alphahat.permu[i] = mean(replicate(m, expr = {
      x = rnorm(group1[i], mu1, sigma1)
      y = rnorm(group2[i], mu2, sigma2)
      z = c(x, y)
      k = sample(1 : (group1[i] + group2[i]), group1[i], replace = FALSE)
      count5test(z[k], z[-k])
    }))
    alphahat.no_permu[i] = mean(replicate(m, expr = {
      x = rnorm(group1[i], mu1, sigma1)
      y = rnorm(group2[i], mu2, sigma2)
      count5test(x, y)
    }))
  }
  data.frame("n1" = group1, "n2" = group2, 
             "perm t1e" = alphahat.permu,
             "no.perm t1e" = alphahat.no_permu)
}



#' @title \code{NN} test statistic calculation
#' @description calculate the statistic of \code{NN} test
#' @import RANN
#' @useDynLib StatComp20084
#' @param z,ix,sizes,k using for statistic calculation
#' @export 
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; 
n2 <- sizes[2]; 
n <- n1 + n2
if(is.vector(z)) 
  z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); 
i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n) 
}





#' @title \code{NN} test p.value calculation
#' @description calculate the p.value of \code{NN} test
#' @import boot
#' @useDynLib StatComp20084
#' @param z,sizes,k using for statistic and p.value calculation
#' @return a list of statistic and p.value
#' @export 
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=999, 
                 sim = "permutation", sizes = sizes,k=k) 
ts <- c(boot.obj$t0,boot.obj$t) 
p.value <- mean(ts>=ts[1]) 
list(statistic=ts[1],p.value=p.value)
}

