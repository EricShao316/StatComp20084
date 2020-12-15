#' @title Ordinary Factor Gradient Descent.
#' @name OFGD
#' @description Factor Gradient Descent considering data X \code{FGD.Ordinary}. Loss Function of Multi-logit model \code{loss.multilogit} and \code{loss.data}. 
#' @description And several functions for simulation setting and accuracy calculation, such as \code{X.generate},\code{trueTheta.generate}.\code{ijS.generate}.
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @useDynLib StatComp20084
#' @examples
#' \dontrun{
#' r1 = 2
#' q1 = 100
#' p1 = 100
#' N1 = 200
#' U1 = matrix(rnorm(p1*2*r1,0,1),p1,2*r1)
#' V1 = matrix(rnorm(q1*2*r1,0,1),q1,2*r1)
#' X1 = X.generate(N1, p1)
#' Theta1 = trueTheta.generate(p1, q1, r1, N1)
#' ijS1 = ijS.generate(N1, p1, q1, X1, Theta1)
#' Theta.hat1 = FGD.Ordinary(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,Data.St = ijS1$St,
#'                           r = 2*r1,tau=1e-5,U = U1,V = V1) $ Ma
#' }
NULL





#' @title Factor Gradient Descent.
#' @description Factor Gradient Descent considering data X. 
#' @param X data matrix (numeric matrix)
#' @param Data.it data it (numeric vector)
#' @param Data.jt data jt (numeric vector)
#' @param Data.St data St (numeric matrix)
#' @param r rank of Theta (integer)
#' @param tau Tolerance (numeric)
#' @param U initail U (numeric matrix)
#' @param V initail V (numeric matrix)
#' @return the estimate of Theta
#' @examples
#' \dontrun{
#' r1 = 2
#' q1 = 100
#' p1 = 100
#' N1 = 200
#' U1 = matrix(rnorm(p1*2*r1,0,1),p1,2*r1)
#' V1 = matrix(rnorm(q1*2*r1,0,1),q1,2*r1)
#' X1 = X.generate(N1, p1)
#' Theta1 = trueTheta.generate(p1, q1, r1, N1)
#' ijS1 = ijS.generate(N1, p1, q1, X1, Theta1)
#' FGD.Ordinary(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,Data.St = ijS1$St,
#'              r = 2*r1,tau=1e-5,U = U1,V = V1)
#' }
#' @export
FGD.Ordinary = function(X = 0, # data matrix
                        Data.it, # data it
                        Data.jt, # data jt
                        Data.St, # data St
                        r, # rank of Theta
                        tau, # Tolerance
                        U, # initial U
                        V # initial V
){
  k = 0 # Times of interation
  K = 50
  f = 200
  f.d = 300
  N = length(Data.it)
  if (is.null(nrow(U))) m = length(U)
  else m = nrow(U)
  if (is.null(nrow(V))) n = length(V)
  else n = nrow(V)
  it = Data.it
  jt = Data.jt
  eit = matrix(0, m, 1)
  ejt = matrix(0, n, 1)
  U.d = matrix(0, m, r)
  V.d = matrix(0, n, r)
  # Regularizing Coefficient
  lambda = sqrt( 10*(m+n)*log(m+n)/(m*n*N) ) / 8 # K = 10
  
  while ((abs(f - f.d) / f.d > tau) & k<30) {
    eta = 1
    f = f.d
    delta.U = - lambda * U
    delta.V = - lambda * V
    
    for (t in 1:N) {
      Q = 0
      weV = matrix(0, n, r)
      weU = matrix(0, m, r)
      J = numeric(0)
      omg = 0
      Wu = matrix(0, m, r)
      Wv = matrix(0, n, r)
      
      eit = matrix(0, m, 1)
      ejt = matrix(0, n, 1)
      eit[it[t],1] = 1
      if (jt[t] != 0) ejt[jt[t],1] = 1
      J = Data.St[t,which(Data.St[t,]!=0)]
      
      for (j in 1:J) {
        ejt.j = matrix(0, n, 1)
        ejt.j[j,1] = 1
        omgu = matrix(exp(- t(U[it[t],]) %*% V[j,] ), m, r)
        omgv = matrix(exp(- t(U[it[t],]) %*% V[j,] ), n, r)
        Wu = Wu + omgu
        Wv = Wv + omgv
        # Sum parts
        weU = weU +  omgu*(eit %*% t(V[j,])) #Drop omg
        weV = weV +  omgv*(ejt.j %*% t(U[it[t],])) #Drop omg
      }
      
      if (jt[t] != 0) delta.U = delta.U - (eit %*% t(V[jt[t],]) - weU / Wu) / N #Drop W
      else delta.U = delta.U - (eit %*% matrix(0, 1, r) - weU / Wu) / N #Drop W
      delta.V = delta.V - (ejt %*% t(U[it[t],]) - weV / Wv) / N #Drop W
    }
    kk = 0
    
    while ((f.d >= f) & (kk < 1e4)) {
      U.d = matrix(0, m, r)
      V.d = matrix(0, n, r)
      U.d = U + eta * delta.U
      V.d = V + eta * delta.V
      
      if (X == 0) LUV.d = loss.multilogit(U.d, V.d, N, it, jt, Data.St)
      else LUV.d = loss.data(X, U.d, V.d, N, it, jt, Data.St, lambda)
      f.d = LUV.d + .5*lambda*sum(U.d^2) + .5*lambda*sum(V.d^2)
      
      # eta.dec = 0.8
      eta = 0.8 * eta
      kk = kk+1
    }
    U = U.d
    V = V.d
    k = k+1
  }
  
  list(Ma = U %*% t(V), f.d = f.d, LUV.d = LUV.d, k = k)
}



#' @title Loss Function with no data matrix.
#' @description calculate the Loss Function with no data matrix in multi-logit model.
#' @param U the estimated U (numeric matrix)
#' @param V the estimated V (numeric matrix)
#' @param N sample size (integer)
#' @param it types (numeric vector)
#' @param jt items (numeric vector)
#' @param St assortment (numeric matrix)
#' @return the value of loss function
#' @examples
#' \dontrun{
#' Loss.value = loss.multilogit(X = X1, U1, V1, N1, ijS1$it, ijS1$jt, ijS1$St)
#' Loss.value
#' }
#' @export
loss.multilogit = function (U, V, N, it, jt, St) {
  Theta = U %*% t(V)
  L = 0
  for (i in 1:N) {
    logit = 1
    for (j in St[i,which(St[i,]!=0)]) {
      logit = logit + exp(Theta[it[i],jt[j]])
    }
    logit = logit * ((jt[i] == 0) + (jt[i] != 0) * exp(Theta[it[i],jt[i]])) ^ (-1)
    L = L + log(logit)
  }
  L = L / N
  return(L)
}



#' @title Loss Function with data matrix.
#' @description calculate the Loss Function with data matrix in multi-logit model.
#' @param X data matrix (numeric matrix)
#' @param U the estimated U (numeric matrix)
#' @param V the estimated V (numeric matrix)
#' @param N sample size (integer)
#' @param it types (numeric vector)
#' @param jt items (numeric vector)
#' @param St assortment (numeric matrix)
#' @param lambda tunning parameter (numeric)
#' @return the value of loss function
#' @examples
#' \dontrun{
#' Loss.value = loss.multilogit(X = X1, U1, V1, N1, ijS1$it, ijS1$jt, ijS1$St, lambda = 0.3)
#' Loss.value
#' }
#' @export
loss.data = function(X, U, V, N, it, jt, St, lambda = 0.3){
  #Theta = U %*% t(V)
  #nuclear = sum(svd(Theta)$d)
  U = X %*% U
  Theta = U %*% t(V)
  L = 0
  for (i in 1:N) {
    logit = 1
    for (j in St[i,]) {
      logit = logit + exp(Theta[i,j])
    }
    
    if (jt[i] != 0) logit = logit / exp(Theta[i,jt[i]])
    else logit = logit
    L = L + log(logit)
  }
  #L = L / N + lambda * nuclear
  L = L / N
  return(L)
}



#' @title Data matrix generated function
#' @description Generate data X.
#' @param N sample size (integer)
#' @param p dimension of features (integer)
#' @return the generated data matrix X
#' @examples
#' \dontrun{
#' X.generate(N = 100, p = 10)
#' }
#' @export
X.generate = function(N, p){
  X = matrix(rnorm(N * p, 0, 1), N, p)
  return(X)
}



#' @title True Theta generated function
#' @description Generate true Theta.
#' @importFrom stats sd
#' @param p number of row (numeric)
#' @param q number of column (numeric)
#' @param r rank of Theta (numeric)
#' @param N sample sizes (numeric)
#' @export
trueTheta.generate = function(p, q, r, N){
  Theta0 = matrix(rnorm(p * q, 0, 1), p, q)
  svd.Theta0 = svd(Theta0)
  len.D = length(svd.Theta0$d)
  Diag = rep(0, len.D)
  Diag[1:r] = svd.Theta0$d[1:r]
  D = diag(Diag)
  U = svd.Theta0$u
  V = svd.Theta0$v
  Theta1 = U %*% D %*% t(V)
  Theta.Star = Theta1 / sd(Theta1)
  for(i in 1:q){
    Theta.Star[sample(1:p,sample(1:floor(p/2))),i] = 0
  }
  return(Theta.Star)
}




#' @title Types, items, assortments generated function
#' @description Generate types, items, assortments.
#' @param N Sample size(time) (numeric)
#' @param p nrow of Theta (numeric)
#' @param q ncol of Theta (numeric)
#' @param X Data Matrix (numeric matrix)
#' @param Theta The generated Theta (numeric matrix)
#' @param K The maximum number of items can be chosen in time t (numeric)
#' @export
ijS.generate = function(N, p, q, X, Theta, K = 50){
  it = sample(1:p, N, replace = TRUE)
  St = matrix(0, N, K) 
  for(i in 1:N){
    # set K = 10
    S.t = sample(1:q, K)
    St[i,1:length(S.t)] = S.t
  }
  pj = matrix(0, N, 1+q)
  jt = numeric(N)
  Theta = X %*% Theta
  for(i in 1:N){
    I = it[i]
    S = St[i,]
    sume = sum(exp(Theta[i,S]))
    pj[i,1] = 1 / sume
    for (j in S) {
      pj[i,j+1] = exp(Theta[i,j]) / sume
    }
    jt[i] = sample(0:q, size = 1, prob = pj[i,])
  }
  list(it = it, jt = jt, St = St)
  
}



#' @title RMSE calculation
#' @param Theta The true Theta
#' @param Theta.hat The estimate of Theta
#' @export
RMSE = function (Theta, Theta.hat) {
  p = nrow(Theta)
  q = ncol(Theta)
  return( sqrt( mean( ( Theta - Theta.hat ) ^ 2 ) ) )
  #return( sqrt( mean( ( Theta - Theta.hat ) ^ 2 ) ) )
}



#' @title FPR, FDR, TPR, TNR, P, R calculation
#' @param Theta.true The true Theta
#' @param Theta.hat The estimate of Theta
#' @export
FPR = function (Theta.true, Theta.hat){
  p = nrow(Theta.true)
  q = ncol(Theta.true)
  Theta.true[which(Theta.true!=0)] = 1
  Theta.hat[which(Theta.hat!=0)] = 1
  TN = length(which(Theta.true+Theta.hat == 0))
  TP = length(which(Theta.true+Theta.hat == 2))
  FN = length(which(Theta.true+Theta.hat == 1 & Theta.true == 1))
  FP = length(which(Theta.true+Theta.hat == 1 & Theta.true == 0))
  list(FPR = FP/(FP+TN), TPR = TP/(TP+FN), P = TP/(TP+FP), R = TP/(TP+FN), FNR = FN/(FN+TP))
  
}