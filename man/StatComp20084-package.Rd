\name{StatComp20084-package}
\alias{StatComp20084-package}
\alias{StatComp20084}
\docType{package}
\title{
   Statitical Computing course project of Shao Jingyu (SA20017032).  
}
\description{
  OFGD and CSFGD algorithms in Assortment problem (calculation and simulation). All the homeworks in the course Statistical Computing.
}
\details{
  All the functions used in the course is now aviliable in this packages.
  Moreover, statitical methods of solving assortment personalization problem, see more in vignettes. 
}
\author{
Shao Jingyu SA20017032

Maintainer: Shao Jingyu <2398799117@qq.com>
}
\references{
  N Kallus, M Udell (2020) Dynamic Assortment Personalization in High Dimensions. Operations Research.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
 r1 = 2
 q1 = 100
 p1 = 100
 N1 = 200
 U1 = matrix(rnorm(p1*2*r1,0,1),p1,2*r1)
 V1 = matrix(rnorm(q1*2*r1,0,1),q1,2*r1)
 X1 = X.generate(N1, p1)
 Theta1 = trueTheta.generate(p1, q1, r1, N1)
 ijS1 = ijS.generate(N1, p1, q1, X1, Theta1)
 FGD.Ordinary(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,
 Data.St = ijS1$St,r = 2*r1,tau=1e-5,U = U1,V = V1) # Ma
 
 FGD.modBeSS(X = X1,Data.it = ijS1$it,Data.jt = ijS1$jt,
 Data.St = ijS1$St, r = 2*r1,U = U1,V = V1) $ Theta.estimate
  }
}
