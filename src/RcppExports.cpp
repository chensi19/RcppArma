// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// grad_desc
arma::mat grad_desc(arma::vec Y, arma::mat X, double step, int dim, int itr, int n1, int n2, arma::vec& Qtrace);
RcppExport SEXP _RcppArma_grad_desc(SEXP YSEXP, SEXP XSEXP, SEXP stepSEXP, SEXP dimSEXP, SEXP itrSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP QtraceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< int >::type itr(itrSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Qtrace(QtraceSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_desc(Y, X, step, dim, itr, n1, n2, Qtrace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppArma_grad_desc", (DL_FUNC) &_RcppArma_grad_desc, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppArma(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
