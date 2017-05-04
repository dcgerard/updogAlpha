#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Functions from prob_functions.cpp ---------------------
Rcpp::NumericVector pbias(Rcpp::NumericVector prob,
                          double bias,
                          double seq_error);

Rcpp::NumericVector dbetabinom_cpp(Rcpp::NumericVector x,
                                   Rcpp::NumericVector size,
                                   double alpha_shape,
                                   double beta_shape,
                                   bool return_log);

Rcpp::NumericVector dbetabinom_mu_rho_cpp(Rcpp::NumericVector x,
                                          Rcpp::NumericVector size, double mu,
                                          double rho, bool return_log);

Rcpp::NumericVector dhyper_cpp(Rcpp::NumericVector x,
                               int m, int n, int k);

arma::Cube<double> get_q_array_cpp(int ploidy);

// Functions from objectives.cpp --------------------------

// Functions from utility.cpp ----------------------------
Rcpp::NumericVector logsumexp(Rcpp::NumericMatrix xx);
