#' @title Competitive Sparse Factor Gradient Descent.
#' @name CSFGD
#' @description Columnwise PDAS considering data X \code{FGD.PDAS}. 
#' @description Competitive Sparse Factor Gradient Descent \code{FGD.modBeSS}, which has high speed in calculation.
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @useDynLib StatComp20084
#' @param p nrow of Theta (numeric)
#' @param X Data matrix (numeric matrix)
#' @param jt item in t (numeric vector)
#' @param Data.St Assortment in t (numeric matrix)
#' @param Theta.j colwise PDAS (numeric)
#' @param U.Ori Initail U (numeric matrix)
#' @param V.Ori Initail V (numeric matrix)
#' @param Theta.1st Prior estimated Theta using to calculate g,h,delta (numeric matrix)
#' @param ActiveSet colwise active set (numeric vector)
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
#' PDAS1 = FGD.PDAS(p = p1, X = X1, jt = ijS1$jt, Data.St = ijS1$St, 
#'                  U.Ori = U1, V.Ori = V1, Theta.1st = U1%*%t(V1), ActiveSet = c(1,3,5,6))
#' PDAS1
#' }
#' @export
FGD.PDAS = function(p, # nrow of Theta
                    X, # Data matrix
                    jt, # item in t
                    Data.St, # Assortment in t
                    Theta.j = 1, # colwise PDAS
                    U.Ori, # Initail U
                    V.Ori, # Initail V
                    Theta.1st, # using to calculate g,h,delta
                    ActiveSet # colwise active set
                    ){
  N = length(jt)
  A1 = rep(1,p)
  A = rep(0,p)
  #ActiveSet = sample(1:p, cardi)
  A[ActiveSet] = A1[ActiveSet] # {1,3,4,6} p=10
  I1 = A1 - A
  ISet = which(I1!=0) 
  # adapt it jt St
  r = ncol(U.Ori)
  U = U.Ori[ActiveSet,] # let col=1
  V = V.Ori[Theta.j,]
  k = 0 # Times of interation
  tau = 1e-5
  f = 200
  f.d = 300
  m = length(ActiveSet)
  n = length(V)
  #it = sample(1:min(N,length(ActiveSet)), N, replace = TRUE)
  it = sample(1:m, N, replace = TRUE)
  jt = sample(0:1, N, replace = TRUE)
  #Data.St = ijS.generate(N, m, n, X[,ActiveSet], U %*% V, K = 1) $ St
  Data.St = matrix(1, N, 1)
  #for (mm in 1:Mmax) {
    # Unit Rank FGD procedure
    eit = matrix(0, m, 1)
    ejt = matrix(0, 1, 1)
    U.d = matrix(0, m, r)
    V.d = matrix(0, n, 1)
    # Regularizing Coefficient
    lambda = sqrt( 10*(m+n)*log(m+n)/(m*n*N) ) / 8 # K = 10
    while ((abs(f - f.d) / f.d > tau) & k<100) {
      k = k + 1
      eta = 1
      f = f.d
      delta.U = - lambda * U
      delta.V = - lambda * V
      
      for (t in 1:N) {
        weV = matrix(0, n, 1)
        weU = matrix(0, m, r)
        J = numeric(0)
        Wu = matrix(0, m, r)
        Wv = matrix(0, n, 1)
        
        eit = matrix(0, m, 1)
        ejt = matrix(0, 1, 1)
        eit[it[t],1] = 1
        if (jt[t] != 0) ejt[jt[t],1] = 1
        #J = Data.St[t,which(Data.St[t,]!=0)]
        J = 1
        
        ejt.j = 1
        #ejt.j[J,1] = 1
        omgu = matrix(exp(- U[it[t],] %*% V ), m, r)
        omgv = matrix(exp(- U[it[t],] %*% V ), n, 1)
        Wu = Wu + omgu
        Wv = Wv + omgv
        # Sum parts
        weU = weU +  omgu * (eit %*% t(V)) #Drop omg
        weV = weV +  omgv * (ejt.j * U[it[t],]) #Drop omg
        
        if (jt[t] != 0) delta.U = delta.U - (eit %*% t(V) - weU / Wu) / N #Drop W
        else delta.U = delta.U - (matrix(0,nrow(weU),ncol(weU)) - weU / Wu) / N #Drop W
        delta.V = delta.V - (ejt * U[it[t],] - weV / Wv) / N #Drop W
      #}
      kk = 0
      
      while ((f.d >= f) & (kk < 1e4)) {
        U.d = matrix(0, m, r)
        V.d = matrix(0, n, 1)
        U.d = U + eta * delta.U
        V.d = V + eta * delta.V
        
        LUV.d = loss.data(X[,ActiveSet], U.d, t(V.d), N, it, jt, Data.St, lambda)
        f.d = LUV.d + .5*lambda*sum(U.d^2) + .5*lambda*sum(V.d^2)
        
        # eta.dec = 0.8
        eta = 0.8 * eta
        kk = kk+1
      }
      U = U.d
      V = V.d
      k = k+1
    }
    Thetaj.estimate = rep(0,p)
    Thetaj.estimate[ActiveSet] = U %*% V
    
    # Now we get the estimate of a col of Theta
    # Then use g,h,delta
    g.colwise = g.calcu(X, Theta = Theta.1st, Theta.j, it, jt, Data.St)
    h.colwise = h.calcu(X, Theta = Theta.1st, Theta.j, it, jt, Data.St)
    delta.colwise = delta.calcu(g = g.colwise, h = h.colwise, Theta = Theta.1st, Theta.j)
  }
  #return(Thetaj.estimate)
  list(Thetaj.estimate = Thetaj.estimate,
       delta.colwise = delta.colwise)
}




#' @title Competitive Sparse Factor Gradient Descent
#' @description Competitive Sparse Factor Gradient Descent considering data X \code{FGD.modBeSS}.
#' @param X data matrix (numeric matrix)
#' @param Data.it data it (numeric vector)
#' @param Data.jt data jt (numeric vector)
#' @param Data.St data St (numeric matrix)
#' @param r rank of Theta (integer)
#' @param U initail U (numeric matrix)
#' @param V initail V (numeric matrix)
#' @param Actint.m m max in BeSS (integer)
#' @param cardi control the level of columnwise sparse (integer)
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
#' Theta.hat1 = FGD.modBeSS(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,Data.St = ijS1$St,
#'                           r = 2*r1,U = U1,V = V1) $ Theta.estimate
#' RMSE(Theta.hat1, trueTheta)
#' }
#' @export
FGD.modBeSS = function(X,
                    Data.it, # data it
                    Data.jt, # data jt
                    Data.St, # data St
                    r, # rank of Theta
                    U, # initial U
                    V, # initial V
                    Actint.m =10, # m max in BeSS
                    cardi = 5
){
  t1 = proc.time()
  # Define
  k = 0 # Times of interation
  kk = 0
  N = nrow(Data.it)
  m = nrow(U)
  n = nrow(V)
  it = Data.it
  jt = Data.jt
  eit = matrix(0, m, 1)
  ejt = matrix(0, n, 1)
  U.d = matrix(0, m, r)
  V.d = matrix(0, n, r)
  cardi = floor(m / 2) # cardinality of active set
  Delta.colwise = matrix(0, m ,n)
  
  ## Theta initailization
  #Theta.estimate = FGD.Ordinary(X,it,jt,Data.St,r,tau=1e-10,U,V)$Ma
  Theta.estimate = U %*% t(V)
  
  ## Active set initailization, K = 10
  A = matrix(0, m, n)
  for (i in 1:n) {
    idx = sample(1:m, sample(2:5))
    A[1:length(idx),i] = idx
  }
  
  
  ## Colwise Update
  for (i in 1:n) {
    pdas = FGD.PDAS(p = m, X = X, jt = Data.jt, 
                    Data.St = Data.St, Theta.j = i, 
                    U.Ori = U, V.Ori = V, Theta.1st = Theta.estimate,
                    ActiveSet = A[which(A[,i]!=0),i])
    Theta.estimate[,i] = pdas $ Thetaj.estimate
    Delta.colwise[,i] = pdas $ delta.colwise
  }
  
  ## Compete Update
  mm = 1
  while (mm <= Actint.m) {
    # Compete on delta to choose one column
    idx = as.vector(vapply(as.data.frame(Delta.colwise),max,numeric(1)))
    idx.max = which(idx==max(idx))
    # Calcu the active set on the choosen col
    kstar = order(Delta.colwise[,idx.max],decreasing = TRUE)[1:cardi]
    A[,idx.max] = rep(0,m)
    A[1:cardi,idx.max] = kstar
    pdas = FGD.PDAS(p = m, X = X, jt = Data.jt, 
                    Data.St = Data.St, Theta.j = idx.max, 
                    U.Ori = U, V.Ori = V, Theta.1st = Theta.estimate,
                    ActiveSet = A[1:cardi, idx.max])
    Theta.estimate[,idx.max] = pdas $ Thetaj.estimate
    Delta.colwise[,idx.max] = pdas $ delta.colwise
    mm = mm + 1
  }
  #return(Theta.estimate)
  t2 = proc.time()
  list(Theta.estimate = Theta.estimate,
       Delta.colwise = Delta.colwise)
}

