// Functions for performing posterior inference in beta-binomial model


#include "updog.h"

//' Posterior inference under the beta-binomial model with outliers. If \code{od_param = 0}, then
//' this reduces to the binomial model with outliers.
//'
//' @param x The observed number of reference alleles.
//' @param n The total number of reads.
//' @param ploidy The ploidy of the species
//' @param prob_geno The distribution of genotypes.
//' @param seq_error The sequencing error rate.
//' @param bias_val The bias parameter. A value of 1 indicates no bias.
//' @param od_param The over-dispersion parameter. Should be in [0, 1). A value of 0 reduces to the
//'     binomial model.
//' @param outlier A logical. Should we include the outlier model (\code{TRUE}) or not (\code{FALSE})?
//' @param out_prop The proportion of observations that are outliers.
//' @param out_mean The mean of outlier distribution.
//' @param out_disp The overdispersion parameter of the outlier distribution.
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector bbpost_double(double x, double n, int ploidy, Rcpp::NumericVector prob_geno,
                                  double seq_error, double bias_val, double od_param,
                                  bool outlier = false, double out_prop = 0.1,
                                  double out_mean = 0.5, double out_disp = 1/3) {
  double tol = 10.0 * DBL_EPSILON;
  // Check input ------------------------------------------------------------------------
  if ((x < 0) | (x > n)) {
    Rcpp::stop("x must be between 0 and n");
  }
  if (ploidy % 2 != 0) {
    Rcpp::stop("ploidy must be even.");
  }
  if (ploidy < 2) {
    Rcpp::stop("ploidy must be at least 2.");
  }
  if (bias_val < 0) {
    Rcpp::stop("bias_val must be greater than 0");
  }
  if ((seq_error < 0) | (seq_error > 1)) {
    Rcpp::stop("seq_error must be between 0 and 1 (inclusive).");
  }
  if ((od_param < 0) | (od_param > 1)) {
    Rcpp::stop("od_param must be between 0 and 1 (inclusive).");
  }
  if ((out_prop < 0) | (out_prop > 1 - tol)) {
    Rcpp::stop("out_prop must be in [0, 1).");
  }
  if ((out_mean < 0) | (out_mean > 1)) {
    Rcpp::stop("out_mean must be between 0 and 1 (inclusive).");
  }
  if ((out_disp < 0) | (out_disp > 1)) {
    Rcpp::stop("out_disp must be between 0 and 1 (inclusive).");
  }

  // get unnormalized probabilities -----------------------------------------------------
  Rcpp::NumericVector pvec = get_pvec(ploidy, bias_val, seq_error);

  Rcpp::NumericVector prob_vec(ploidy + 1);
  for (int i = 0; i < ploidy + 1; i++) {
    if (prob_geno(i) > tol) {
      prob_vec(i) = dbetabinom_mu_rho_cpp_double(x, n, pvec(i), od_param, true) +
        log(prob_geno(i));
    } else {
      prob_vec(i) = R_NegInf;
    }
  }

  // outlier model if outlier = true by point-wise log-sum-exp --------------------------
  if (outlier & (out_prop > tol)) {
    double out_dense = std::log(out_prop) + dbetabinom_mu_rho_cpp_double(x, n, out_mean, out_disp, true);
    Rcpp::NumericVector outlier_vec(ploidy + 1, out_dense);
    prob_vec = prob_vec + log(1 - out_prop);
    Rcpp::NumericVector max_el = Rcpp::pmax(prob_vec, outlier_vec);
    prob_vec = Rcpp::log(Rcpp::exp(prob_vec - max_el) + Rcpp::exp(outlier_vec - max_el)) + max_el;
  }

  // log-sum-exp probs ------------------------------------------------------------------
  double lmax = Rcpp::max(prob_vec);
  double lsum = std::log(Rcpp::sum(Rcpp::exp(prob_vec - lmax))) + lmax;
  prob_vec = Rcpp::exp(prob_vec - lsum);
  return prob_vec;
}

//' Posterior inference for each individual.
//'
//' @inheritParams bbpost_double
//' @param ocounts The vector of counts of the reference allele of the offspring.
//' @param osize The vector of total reads of the offspring.
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix bbpost_tot(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                               int ploidy, Rcpp::NumericVector prob_geno,
                               double seq_error, double bias_val, double od_param,
                               bool outlier = false, double out_prop = 0.1,
                               double out_mean = 0.5, double out_disp = 1/3) {
  if (ocounts.size() != osize.size()) {
    Rcpp::stop("ocounts and osize need to be of the same length");
  }

  Rcpp::NumericMatrix prob_mat(ocounts.size(), ploidy + 1);
  for (int i = 0; i < ocounts.size(); i++) {
    prob_mat(i, Rcpp::_) = bbpost_double(ocounts(i), osize(i), ploidy, prob_geno,
             seq_error, bias_val, od_param,
             outlier, out_prop, out_mean, out_disp);
  }
  return(prob_mat);
}



