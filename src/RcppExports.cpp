// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// LCS
int LCS(std::string x, std::string y);
RcppExport SEXP _sigminer_helper_LCS(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(LCS(x, y));
    return rcpp_result_gen;
END_RCPP
}
// LCSMatrix
IntegerMatrix LCSMatrix(StringVector x, StringVector y, bool match);
RcppExport SEXP _sigminer_helper_LCSMatrix(SEXP xSEXP, SEXP ySEXP, SEXP matchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< StringVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type match(matchSEXP);
    rcpp_result_gen = Rcpp::wrap(LCSMatrix(x, y, match));
    return rcpp_result_gen;
END_RCPP
}
// pairScoreVector
int pairScoreVector(NumericVector x, NumericVector y, int x_max, int y_max);
RcppExport SEXP _sigminer_helper_pairScoreVector(SEXP xSEXP, SEXP ySEXP, SEXP x_maxSEXP, SEXP y_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< int >::type y_max(y_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(pairScoreVector(x, y, x_max, y_max));
    return rcpp_result_gen;
END_RCPP
}
// pairScoreMatrix
NumericMatrix pairScoreMatrix(NumericMatrix x, NumericMatrix y, int x_max, int y_max);
RcppExport SEXP _sigminer_helper_pairScoreMatrix(SEXP xSEXP, SEXP ySEXP, SEXP x_maxSEXP, SEXP y_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< int >::type y_max(y_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(pairScoreMatrix(x, y, x_max, y_max));
    return rcpp_result_gen;
END_RCPP
}
// pairScoreSimpleVector
int pairScoreSimpleVector(NumericVector x, NumericVector y, int max);
RcppExport SEXP _sigminer_helper_pairScoreSimpleVector(SEXP xSEXP, SEXP ySEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(pairScoreSimpleVector(x, y, max));
    return rcpp_result_gen;
END_RCPP
}
// pairScoreSimpleMatrix
NumericMatrix pairScoreSimpleMatrix(NumericMatrix x, NumericMatrix y, int max);
RcppExport SEXP _sigminer_helper_pairScoreSimpleMatrix(SEXP xSEXP, SEXP ySEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(pairScoreSimpleMatrix(x, y, max));
    return rcpp_result_gen;
END_RCPP
}
// getScoreMatrix
IntegerMatrix getScoreMatrix(IntegerMatrix indexMat, IntegerMatrix subMat, int bSize, bool like, bool verbose);
RcppExport SEXP _sigminer_helper_getScoreMatrix(SEXP indexMatSEXP, SEXP subMatSEXP, SEXP bSizeSEXP, SEXP likeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type indexMat(indexMatSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type subMat(subMatSEXP);
    Rcpp::traits::input_parameter< int >::type bSize(bSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type like(likeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(getScoreMatrix(indexMat, subMat, bSize, like, verbose));
    return rcpp_result_gen;
END_RCPP
}
// getScoreMatrixRect
IntegerMatrix getScoreMatrixRect(IntegerMatrix indexMat1, IntegerMatrix indexMat2, IntegerMatrix subMat, bool like, bool verbose);
RcppExport SEXP _sigminer_helper_getScoreMatrixRect(SEXP indexMat1SEXP, SEXP indexMat2SEXP, SEXP subMatSEXP, SEXP likeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type indexMat1(indexMat1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indexMat2(indexMat2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type subMat(subMatSEXP);
    Rcpp::traits::input_parameter< bool >::type like(likeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(getScoreMatrixRect(indexMat1, indexMat2, subMat, like, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sigminer_helper_LCS", (DL_FUNC) &_sigminer_helper_LCS, 2},
    {"_sigminer_helper_LCSMatrix", (DL_FUNC) &_sigminer_helper_LCSMatrix, 3},
    {"_sigminer_helper_pairScoreVector", (DL_FUNC) &_sigminer_helper_pairScoreVector, 4},
    {"_sigminer_helper_pairScoreMatrix", (DL_FUNC) &_sigminer_helper_pairScoreMatrix, 4},
    {"_sigminer_helper_pairScoreSimpleVector", (DL_FUNC) &_sigminer_helper_pairScoreSimpleVector, 3},
    {"_sigminer_helper_pairScoreSimpleMatrix", (DL_FUNC) &_sigminer_helper_pairScoreSimpleMatrix, 3},
    {"_sigminer_helper_getScoreMatrix", (DL_FUNC) &_sigminer_helper_getScoreMatrix, 5},
    {"_sigminer_helper_getScoreMatrixRect", (DL_FUNC) &_sigminer_helper_getScoreMatrixRect, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_sigminer_helper(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
