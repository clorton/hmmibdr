// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// hmmidm_cpp
int hmmidm_cpp(Rcpp::List param_list);
RcppExport SEXP _hmmidm_hmmidm_cpp(SEXP param_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param_list(param_listSEXP);
    rcpp_result_gen = Rcpp::wrap(hmmidm_cpp(param_list));
    return rcpp_result_gen;
END_RCPP
}
// IdmRcppTest
DataFrame IdmRcppTest(const DataFrame& dfin);
RcppExport SEXP _hmmidm_IdmRcppTest(SEXP dfinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type dfin(dfinSEXP);
    rcpp_result_gen = Rcpp::wrap(IdmRcppTest(dfin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hmmidm_hmmidm_cpp", (DL_FUNC) &_hmmidm_hmmidm_cpp, 1},
    {"_hmmidm_IdmRcppTest", (DL_FUNC) &_hmmidm_IdmRcppTest, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_hmmidm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
