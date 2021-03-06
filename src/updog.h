#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Functions from prob_functions.cpp ---------------------
Rcpp::NumericVector pbias(Rcpp::NumericVector prob,
                          double bias,
                          double seq_error);

double pbias_double(double prob, double bias, double seq_error);

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
Rcpp::NumericVector get_prob_geno(int ploidy, std::string model, int p1geno, int p2geno, double allele_freq);

double dbetabinom_mu_rho_cpp_double(double x, double size, double mu,
                                    double rho, bool return_log);

// Functions from objectives.cpp --------------------------
Rcpp::NumericVector obj_offspring_vec(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                                      int ploidy, Rcpp::NumericVector prob_geno,
                                      double bias_val, double seq_error,
                                      double od_param,
                                      bool outlier, double out_prop,
                                      double out_mean, double out_disp);

// Functions from utility.cpp ----------------------------
Rcpp::NumericVector logsumexp(Rcpp::NumericMatrix xx);
Rcpp::NumericVector get_pvec(int ploidy, double bias_val, double seq_error);
Rcpp::NumericVector colSums_cpp(Rcpp::NumericMatrix x);
double expit(double x);

// Functions from gradients.cpp --------------------------
double dbeta_dprop(double x, double n, double xi, double tau);
double dbeta_dh(double x, double n, double xi, double h);
double dh_dtau(double tau);
double dbeta_deps(double x, double n, double d, double eps, double p, double tau);
double dbeta_dd(double x, double n, double d, double eps, double p, double tau);
double dbeta_dtau(double x, double n, double d, double eps, double p, double tau);
double dbeta_ds(double x, double n, double s, double ell, double p, double h);
double dbeta_dl(double x, double n, double d, double ell, double p, double h);
double dbeta_dr_ell(double x, double n, double d, double ell, double p, double r);
