#' @title Sparse Factor Gradient Descent.
#' @name SFGD
#' @description Factor Gradient Descent considering data X \code{FGD.BeSS}.
#' @description And several functions for 1-2 order gradient calculation, \code{g.calcu},\code{h.calcu}.\code{delta.calcu}.
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
#' Theta.hat1 = FGD.BeSS(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,Data.St = ijS1$St,
#'                           r = 2*r1,tau=1e-5,U = U1,V = V1) $ Theta.estimate
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
#' @param Actint.m m max in BeSS (integer)
#' @return the estimate of Theta and error.
#' @examples
#' \dontrun{
#' FGD.BeSS(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,Data.St = ijS1$St,
#'              r = 2*r1,tau=1e-5,U = U1,V = V1)
#' }
#' @export
FGD.BeSS = function(X,
                    Data.it, # data it
                    Data.jt, # data jt
                    Data.St, # data St
                    r, # rank of Theta
                    tau, # Tolerance
                    U, # initial U
                    V, # initial V
                    Actint.m =100 # m max in BeSS
){
  t1 = proc.time()
  # Define
  k = 0 # Times of interation
  kk = 0
  f = 200
  f.d = 300
  N = nrow(Data.it)
  m = nrow(U)
  n = nrow(V)
  it = Data.it
  jt = Data.jt
  eit = matrix(0, m, 1)
  ejt = matrix(0, n, 1)
  U.d = matrix(0, m, r)
  V.d = matrix(0, n, r)
  cardi = m - 2 # cardinality of active set
  
  ## Theta initailization
  Theta.estimate = FGD.Ordinary(X,it,jt,Data.St,r,tau=1e-10,U,V)$Ma
  
  ## Active set initailization
  A = matrix(0, m, n)
  errA = matrix(1, m, n)
  for (i in 1:n) {
    idx = sample(1:m, sample(1:(m-2)))
    A[1:length(idx),i] = idx
  }
  
  ## BeSS-FGD
  mm = 0
  while ((mm <= Actint.m) & (sum(abs(errA)) != 0)) {
    # estimate Theta by column, run FGD on active set
    k = k + 1
    g.colwise = matrix(0, m, n)
    h.colwise = matrix(0, m, n)
    delta.colwise = matrix(0, m, n)
    delta = matrix(0, m, n)
    
    for (i in 1:n) {
      idx = A[which(A[,i]!=0),i]
      if (is.null(nrow(U[idx,])))
        U.i = t(as.matrix(U[idx,]))
      else
        U.i = U[idx,]
      it = sample(1:length(idx), N, replace = TRUE)
      X.i = as.matrix(X[,idx])
      Theta.estimate[-idx,i] = 0
      Theta.estimate[idx,i] = FGD.Ordinary(X.i, it, jt,
                                  Data.St,r,tau=1e-10,U.i,V)$Ma[,i]
      g.colwise[,i] = g.calcu(X, Theta = Theta.estimate, Theta.j = i, it, jt, Data.St)
      h.colwise[,i] = h.calcu(X, Theta = Theta.estimate, Theta.j = i, it, jt, Data.St)
      delta.colwise[,i] = delta.calcu(g = g.colwise[,i], h = h.colwise[,i], Theta = Theta.estimate, Theta.j = i)
      delta[,i] = sort(delta.colwise[,i], decreasing = TRUE)
    }
    
    # Type = PDAS
    index = matrix(0, m, n)
    
    # Update Active Set
    for(i in 1:n){
      windex = which(delta.colwise[,i] >= delta[cardi,i]) 
      index[1:length(windex),i] = windex
    }
    
    errA = A - index
    A = index
    mm = mm + 1
  }
  t2 = proc.time()
  
  list(Theta.estimate = Theta.estimate, 
       process.time = t2 - t1,
       k = k - 1,
       errA = errA)
}



#' @title in.St
#' @description in.St used to calculate g and h 
#' @param X,ThetaCol,StRow,itRow,St A preperation for calculating g and h.
in.St = function(X, ThetaCol, StRow, itRow, St) {
  n = nrow(St)
  s = 1
  for (j in StRow){
    s = s + exp(t(X[itRow,]) %*% ThetaCol)
  }
  return(s)
}

#' @title Gradient Calculation
#' @description Gradient Calculation 
#' @param X,Theta,Theta.j,it,jt,St Are using to calculate Gradient.
g.calcu = function(X, Theta, Theta.j, it, jt, St) {
  N = length(it)
  m = ncol(X)
  g = matrix(0,m,1)
  for(i in 1:N){
    inst = in.St(X, ThetaCol = Theta[,Theta.j], StRow = St[i,], itRow = it[i], St = St)
    if (jt[i]!=0) {
      g = g - exp(t(X[it[i],]) %*% (Theta[,Theta.j] + 2 * Theta[,jt[i]])) * X[i,] / (inst)^3
    }
    else{
      g = g - exp(t(X[it[i],]) %*% (Theta[,Theta.j] + 2 * matrix(0, m, 1))) * X[i,] / (inst)^3
    }
    #inst = in.St(X, ThetaCol = Theta[,Theta.j], StRow = St[i,], itRow = it[i], St = St)
    #g = g - exp(t(X[it[i],]) %*% (Theta[,Theta.j] + 2 * Theta[,jt[i]])) * X[i,] / (inst)^3
  }
  return(g / N)
}


#' @title Hessian Calculation
#' @description Hessian Calculation
#' @param X,Theta,Theta.j,it,jt,St Are using to calculate diag Hessian.
h.calcu = function(X, Theta, Theta.j, it, jt, St) {
  N = length(it)
  m = ncol(X)
  h = matrix(0,m,1)
  for(i in 1:N){
    inst = in.St(X, ThetaCol = Theta[,Theta.j], StRow = St[i,], itRow = it[i], St = St)
    if (jt[i]!=0){
      expst = exp(t(X[it[i],]) %*% (Theta[,Theta.j] + 2 * Theta[,jt[i]]))
    }
    else expst = exp(t(X[it[i],]) %*% (Theta[,Theta.j] + 2 * matrix(0, m, 1)))
    #expst = exp(t(X[it[i],]) %*% (Theta[,Theta.j] + 2 * Theta[,jt[i]]))
    h = h + (-expst * X[it[i],] + 3 * expst * inst^2 * exp( t(X[it[i],]) %*% Theta[,Theta.j] ) * X[i,]) / (inst)^6
  }
  return(h)
}


#' @title Delta Calculation
#' @description Delta Calculation
#' @param g,h,Theta,Theta.j Are using to calculate Delta.
#' @examples 
#' \dontrun{
#' g1 = g.calcu(X1,Theta1,Theta.j=5,ijS1$it,ijS1$jt,ijS1$St)
#' h1 = h.calcu(X1,Theta1,Theta.j=5,ijS1$it,ijS1$jt,ijS1$St)
#' delta1 = delta.calcu(g1,h1,Theta1,Theta.j=5)
#' testgh1 = data.frame(g1 = g1, h1 = h1, delta1 = delta1)
#' testgh1
#' }
delta.calcu = function(g, h, Theta, Theta.j) {
  gamma.j = - g / h
  delta.j = .5 * h * (Theta[,Theta.j] + gamma.j) ^ 2
  return(delta.j)
}

