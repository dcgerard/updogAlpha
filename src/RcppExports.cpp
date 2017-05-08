// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dbeta_dprop
double dbeta_dprop(double x, double n, double xi, double tau);
RcppExport SEXP updog_dbeta_dprop(SEXP xSEXP, SEXP nSEXP, SEXP xiSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_dprop(x, n, xi, tau));
    return rcpp_result_gen;
END_RCPP
}
// dbeta_dh
double dbeta_dh(double x, double n, double xi, double h);
RcppExport SEXP updog_dbeta_dh(SEXP xSEXP, SEXP nSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_dh(x, n, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// dxi_df
double dxi_df(double d, double f);
RcppExport SEXP updog_dxi_df(SEXP dSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(dxi_df(d, f));
    return rcpp_result_gen;
END_RCPP
}
// df_deps
double df_deps(double eps, double p);
RcppExport SEXP updog_df_deps(SEXP epsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(df_deps(eps, p));
    return rcpp_result_gen;
END_RCPP
}
// deps_dell
double deps_dell(double ell);
RcppExport SEXP updog_deps_dell(SEXP ellSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type ell(ellSEXP);
    rcpp_result_gen = Rcpp::wrap(deps_dell(ell));
    return rcpp_result_gen;
END_RCPP
}
// expit
double expit(double x);
RcppExport SEXP updog_expit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(expit(x));
    return rcpp_result_gen;
END_RCPP
}
// dbeta_dl
double dbeta_dl(double x, double n, double d, double ell, double p, double h);
RcppExport SEXP updog_dbeta_dl(SEXP xSEXP, SEXP nSEXP, SEXP dSEXP, SEXP ellSEXP, SEXP pSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_dl(x, n, d, ell, p, h));
    return rcpp_result_gen;
END_RCPP
}
// dxi_dd
double dxi_dd(double d, double f);
RcppExport SEXP updog_dxi_dd(SEXP dSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(dxi_dd(d, f));
    return rcpp_result_gen;
END_RCPP
}
// dbeta_dd
double dbeta_dd(double x, double n, double d, double ell, double p, double h);
RcppExport SEXP updog_dbeta_dd(SEXP xSEXP, SEXP nSEXP, SEXP dSEXP, SEXP ellSEXP, SEXP pSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_dd(x, n, d, ell, p, h));
    return rcpp_result_gen;
END_RCPP
}
// dbeta_dh_ell
double dbeta_dh_ell(double x, double n, double d, double ell, double p, double h);
RcppExport SEXP updog_dbeta_dh_ell(SEXP xSEXP, SEXP nSEXP, SEXP dSEXP, SEXP ellSEXP, SEXP pSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_dh_ell(x, n, d, ell, p, h));
    return rcpp_result_gen;
END_RCPP
}
// grad_offspring
Rcpp::NumericVector grad_offspring(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize, int ploidy, int p1geno, int p2geno, double bias_val, double seq_error, double od_param, bool outlier, double out_prop, double out_mean, double out_disp);
RcppExport SEXP updog_grad_offspring(SEXP ocountsSEXP, SEXP osizeSEXP, SEXP ploidySEXP, SEXP p1genoSEXP, SEXP p2genoSEXP, SEXP bias_valSEXP, SEXP seq_errorSEXP, SEXP od_paramSEXP, SEXP outlierSEXP, SEXP out_propSEXP, SEXP out_meanSEXP, SEXP out_dispSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ocounts(ocountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type osize(osizeSEXP);
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
    rcpp_result_gen = Rcpp::wrap(grad_offspring(ocounts, osize, ploidy, p1geno, p2geno, bias_val, seq_error, od_param, outlier, out_prop, out_mean, out_disp));
    return rcpp_result_gen;
END_RCPP
}
// obj_offspring_vec
Rcpp::NumericVector obj_offspring_vec(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize, int ploidy, int p1geno, int p2geno, double bias_val, double seq_error, double od_param, bool outlier, double out_prop, double out_mean, double out_disp);
RcppExport SEXP updog_obj_offspring_vec(SEXP ocountsSEXP, SEXP osizeSEXP, SEXP ploidySEXP, SEXP p1genoSEXP, SEXP p2genoSEXP, SEXP bias_valSEXP, SEXP seq_errorSEXP, SEXP od_paramSEXP, SEXP outlierSEXP, SEXP out_propSEXP, SEXP out_meanSEXP, SEXP out_dispSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ocounts(ocountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type osize(osizeSEXP);
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
    rcpp_result_gen = Rcpp::wrap(obj_offspring_vec(ocounts, osize, ploidy, p1geno, p2geno, bias_val, seq_error, od_param, outlier, out_prop, out_mean, out_disp));
    return rcpp_result_gen;
END_RCPP
}
// obj_offspring
double obj_offspring(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize, int ploidy, int p1geno, int p2geno, double bias_val, double seq_error, double od_param, bool outlier, double out_prop, double out_mean, double out_disp);
RcppExport SEXP updog_obj_offspring(SEXP ocountsSEXP, SEXP osizeSEXP, SEXP ploidySEXP, SEXP p1genoSEXP, SEXP p2genoSEXP, SEXP bias_valSEXP, SEXP seq_errorSEXP, SEXP od_paramSEXP, SEXP outlierSEXP, SEXP out_propSEXP, SEXP out_meanSEXP, SEXP out_dispSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ocounts(ocountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type osize(osizeSEXP);
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
// obj_offspring_weights
double obj_offspring_weights(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize, Rcpp::NumericVector weight_vec, int ploidy, int p1geno, int p2geno, double bias_val, double seq_error, double od_param, bool outlier, double out_prop, double out_mean, double out_disp);
RcppExport SEXP updog_obj_offspring_weights(SEXP ocountsSEXP, SEXP osizeSEXP, SEXP weight_vecSEXP, SEXP ploidySEXP, SEXP p1genoSEXP, SEXP p2genoSEXP, SEXP bias_valSEXP, SEXP seq_errorSEXP, SEXP od_paramSEXP, SEXP outlierSEXP, SEXP out_propSEXP, SEXP out_meanSEXP, SEXP out_dispSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ocounts(ocountsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type osize(osizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weight_vec(weight_vecSEXP);
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
    rcpp_result_gen = Rcpp::wrap(obj_offspring_weights(ocounts, osize, weight_vec, ploidy, p1geno, p2geno, bias_val, seq_error, od_param, outlier, out_prop, out_mean, out_disp));
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
// pbias_double
double pbias_double(double prob, double bias, double seq_error);
RcppExport SEXP updog_pbias_double(SEXP probSEXP, SEXP biasSEXP, SEXP seq_errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type bias(biasSEXP);
    Rcpp::traits::input_parameter< double >::type seq_error(seq_errorSEXP);
    rcpp_result_gen = Rcpp::wrap(pbias_double(prob, bias, seq_error));
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
// dbetabinom_mu_rho_cpp_double
double dbetabinom_mu_rho_cpp_double(double x, double size, double mu, double rho, bool return_log);
RcppExport SEXP updog_dbetabinom_mu_rho_cpp_double(SEXP xSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP return_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type return_log(return_logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom_mu_rho_cpp_double(x, size, mu, rho, return_log));
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
// get_pvec
Rcpp::NumericVector get_pvec(int ploidy, double bias_val, double seq_error);
RcppExport SEXP updog_get_pvec(SEXP ploidySEXP, SEXP bias_valSEXP, SEXP seq_errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< double >::type bias_val(bias_valSEXP);
    Rcpp::traits::input_parameter< double >::type seq_error(seq_errorSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pvec(ploidy, bias_val, seq_error));
    return rcpp_result_gen;
END_RCPP
}
