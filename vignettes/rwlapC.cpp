#include <Rcpp.h>
#include <iostream>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
double dlapC(double x){
	double y=0;
	return(0.5 * exp(-abs(x)));
}



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
