// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dlapC
double dlapC(double x);
RcppExport SEXP _StatComp20084_dlapC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(dlapC(x));
    return rcpp_result_gen;
END_RCPP
}
// rwlapC
List rwlapC(int N, double sigma);
RcppExport SEXP _StatComp20084_rwlapC(SEXP NSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rwlapC(N, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp20084_dlapC", (DL_FUNC) &_StatComp20084_dlapC, 1},
    {"_StatComp20084_rwlapC", (DL_FUNC) &_StatComp20084_rwlapC, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp20084(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
