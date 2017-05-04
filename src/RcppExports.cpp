// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// obj_offspring
double obj_offspring(Rcpp::IntegerVector ocounts, Rcpp::IntegerVector osize, int ploidy, int p1geno, int p2geno, double bias_val, double seq_error, double od_param, bool outlier, double out_prop, double out_mean, double out_disp);
RcppExport SEXP updog_obj_offspring(SEXP ocountsSEXP, SEXP osizeSEXP, SEXP ploidySEXP, SEXP p1genoSEXP, SEXP p2genoSEXP, SEXP bias_valSEXP, SEXP seq_errorSEXP, SEXP od_paramSEXP, SEXP outlierSEXP, SEXP out_propSEXP, SEXP out_meanSEXP, SEXP out_dispSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ocounts(ocountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type osize(osizeSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< int >::type p1geno(p1genoSEXP);
    Rcpp::traits::input_parameter< int >::type p2geno(p2genoSEXP);
    Rcpp::traits::input_parameter< double >::type bias_val(bias_valSEXP);
    Rcpp::traits::input_parameter< double >::type seq_error(seq_errorSEXP);
    Rcpp::traits::input_parameter< double >::type od_param(od_paramSEXP);
    Rcpp::traits::input_parameter< bool >::type outlier(outlierSEXP);
    Rcpp::traits::input_parameter< double >::type out_prop(out_propSEXP);
    Rcpp::traits::input_parameter< double >::type out_mean(out_meanSEXP);
    Rcpp::traits::input_parameter< double >::type out_disp(out_dispSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_offspring(ocounts, osize, ploidy, p1geno, p2geno, bias_val, seq_error, od_param, outlier, out_prop, out_mean, out_disp));
    return rcpp_result_gen;
END_RCPP
}
// pbias
Rcpp::NumericVector pbias(Rcpp::NumericVector prob, double bias, double seq_error);
RcppExport SEXP updog_pbias(SEXP probSEXP, SEXP biasSEXP, SEXP seq_errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type bias(biasSEXP);
    Rcpp::traits::input_parameter< double >::type seq_error(seq_errorSEXP);
    rcpp_result_gen = Rcpp::wrap(pbias(prob, bias, seq_error));
    return rcpp_result_gen;
END_RCPP
}
// dbetabinom_cpp
Rcpp::NumericVector dbetabinom_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, double alpha_shape, double beta_shape, bool return_log);
RcppExport SEXP updog_dbetabinom_cpp(SEXP xSEXP, SEXP sizeSEXP, SEXP alpha_shapeSEXP, SEXP beta_shapeSEXP, SEXP return_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_shape(alpha_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type beta_shape(beta_shapeSEXP);
    Rcpp::traits::input_parameter< bool >::type return_log(return_logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom_cpp(x, size, alpha_shape, beta_shape, return_log));
    return rcpp_result_gen;
END_RCPP
}
// dbetabinom_mu_rho_cpp
Rcpp::NumericVector dbetabinom_mu_rho_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, double mu, double rho, bool return_log);
RcppExport SEXP updog_dbetabinom_mu_rho_cpp(SEXP xSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP return_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type return_log(return_logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom_mu_rho_cpp(x, size, mu, rho, return_log));
    return rcpp_result_gen;
END_RCPP
}
// dhyper_cpp
Rcpp::NumericVector dhyper_cpp(Rcpp::NumericVector x, int m, int n, int k);
RcppExport SEXP updog_dhyper_cpp(SEXP xSEXP, SEXP mSEXP, SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(dhyper_cpp(x, m, n, k));
    return rcpp_result_gen;
END_RCPP
}
// get_q_array_cpp
arma::Cube<double> get_q_array_cpp(int ploidy);
RcppExport SEXP updog_get_q_array_cpp(SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(get_q_array_cpp(ploidy));
    return rcpp_result_gen;
END_RCPP
}
// logsumexp
Rcpp::NumericVector logsumexp(Rcpp::NumericMatrix xx);
RcppExport SEXP updog_logsumexp(SEXP xxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xx(xxSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp(xx));
    return rcpp_result_gen;
END_RCPP
}
