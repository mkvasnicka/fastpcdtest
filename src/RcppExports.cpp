// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fastpcdtest_balanced_
double fastpcdtest_balanced_(NumericMatrix M);
RcppExport SEXP fastpcdtest_fastpcdtest_balanced_(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    __result = Rcpp::wrap(fastpcdtest_balanced_(M));
    return __result;
END_RCPP
}
// fastpcdtest_unbalanced_
NumericVector fastpcdtest_unbalanced_(NumericMatrix M, int min_common);
RcppExport SEXP fastpcdtest_fastpcdtest_unbalanced_(SEXP MSEXP, SEXP min_commonSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type min_common(min_commonSEXP);
    __result = Rcpp::wrap(fastpcdtest_unbalanced_(M, min_common));
    return __result;
END_RCPP
}