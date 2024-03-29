% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FGD_Ordinary.R
\name{loss.data}
\alias{loss.data}
\title{Loss Function with data matrix.}
\usage{
loss.data(X, U, V, N, it, jt, St, lambda = 0.3)
}
\arguments{
\item{X}{data matrix (numeric matrix)}

\item{U}{the estimated U (numeric matrix)}

\item{V}{the estimated V (numeric matrix)}

\item{N}{sample size (integer)}

\item{it}{types (numeric vector)}

\item{jt}{items (numeric vector)}

\item{St}{assortment (numeric matrix)}

\item{lambda}{tunning parameter (numeric)}
}
\value{
the value of loss function
}
\description{
calculate the Loss Function with data matrix in multi-logit model.
}
\examples{
\dontrun{
Loss.value = loss.multilogit(X = X1, U1, V1, N1, ijS1$it, ijS1$jt, ijS1$St, lambda = 0.3)
Loss.value
}
}
