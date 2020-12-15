#include <Rcpp.h>
using namespace Rcpp;

//' @title Laplace density function using Rcpp
//' @description Laplace density function using Rcpp
//' @param x the sample points
//' @return the value of density function in x
//' @examples
//' \dontrun{
//' dlap <- dlapC(5)
//' dlap
//' }
//' @export
// [[Rcpp::export]]
double dlapC(double x){
	double y=0;
	return(0.5 * exp(-abs(x)));
}

#include <Rcpp.h>
using namespace Rcpp;
//' @title A Random Walk sampler using Rcpp
//' @description A Random Walk sampler using Rcpp
//' @param N the length of chains
//' @param sigma the variance of normal distribution
//' @return the generated chain
//' @examples
//' \dontrun{
//' rwlap <- rwlapC(5000, 2)
//' plot(1:5000, rwlap$x, ylab = "x", 
//' main = expression(paste(sigma, '= 2')), type = "l")
//' }
//' @export
// [[Rcpp::export]]
List rwlapC(int N, double sigma){
	NumericVector x(N);
	NumericVector u(N);
	double y = 0;
	int k = 0;
	List lapC(0);
	x[0] = 0;
	u = runif(N);
	for(int i = 1; i < N; i++){
		y = rnorm(1, x[i-1], sigma)[0];
		if(u[i] <= dlapC(y) / dlapC(x[i-1])){
			x[i] = y;
		}
		else {
			x[i] = x[i-1];
			k += 1;
		}
	}
	lapC["x"] = x;
	lapC["k"] = k;
	return(lapC);
}
