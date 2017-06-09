////////
// This file contains updog objective functions
////////

#include "updog.h"

//' Vector of objective functions for offspring.
//'
//' @param ocounts The observed counts of the refernce
//'     allele for each individual.
//' @param osize The observed number of reads for each
//'     individuals.
//' @param ploidy An integer. The ploidy of the species. This is assumed
//'     to be the same for each individual.
//' @param prob_geno The allele frequencies of the genotypes. See \code{\link{get_prob_geno}}.
//' @param bias_val The bias parameter. A value of 1 means there is no bias
//'     toward one allele or the other. A value less than one indicates a bias
//'     toward the reference allele. A value greater than one indicates a bias
//'     toward the non-reference allele.
//' @param seq_error The sequencing error rate.
//' @param od_param The overdispersion parameter in the beta-binomial model
//'     for the OK counts. When this is zero, this resorts to the binomial
//'     model for counts.
//' @param outlier A logical. Should we include an outlier model (\code{TRUE})
//'     or not (\code{FALSE})? Defaults to \code{FALSE}.
//' @param out_prop The proportion of points that are outliers. Defaults
//'     (quite arbitrarily) to \code{0.01}.
//' @param out_mean The mean of beta-binomial for the outlier distribution.
//'     Defaults to \code{0.5}.
//' @param out_disp The overdispersion parameter of the outlier distribution.
//'     Defaults to \code{1/3}, which corresponds to a uniform distribution
//'     for the underlying beta when \code{out_mean = 0.5}.
//'
//' @author David Gerard
//'
//'
//' @seealso \code{\link{up_bb_obj}}.
//'
// [[Rcpp::export]]
Rcpp::NumericVector obj_offspring_vec(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                                      int ploidy, Rcpp::NumericVector prob_geno,
                                      double bias_val = 1, double seq_error = 0,
                                      double od_param = 0,
                                      bool outlier = false, double out_prop = 0.01,
                                      double out_mean = 0.5, double out_disp = 1.0 / 3.0) {
  double tol = 2 * DBL_EPSILON;

  // check input ---------------------------------------------------
  if (ocounts.size() != osize.size()) {
    Rcpp::stop("ocounts and osize must have the same length.");
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

  // calculate the possible means proportions based on 0:ploidy / ploidy,
  // seq_error, and bias_val.
  Rcpp::NumericVector prob = get_pvec(ploidy, bias_val, seq_error);

  // Calculate probabilities for OK points --------------------------
  // only need last column if outlier = true ------------------------

  // count how many non-zero values in prob_geno values there are
  int colnum = 0;
  for (int i = 0; i < ploidy + 1; i++) {
    if (prob_geno(i) > tol) {
      colnum++;
    }
  }
  // add one for outlier
  if (outlier) {
    colnum++;
  }

  Rcpp::NumericMatrix logprobs(ocounts.size(), colnum);

  // calcluate log mixture component probabilities -----------------
  int colindex = 0;
  for (int i = 0; i < ploidy + 1; i++) {
    if (prob_geno(i) > tol) {
      Rcpp::NumericMatrix::Column zzcol = logprobs(Rcpp::_, colindex); // reference to colindexth column
      zzcol = dbetabinom_mu_rho_cpp(ocounts, osize, prob(i), od_param, true) +
	log(prob_geno(i));
      if (outlier) {
        zzcol = zzcol + log(1 - out_prop);
      }
      colindex++;
    }
  }

  // Another beta binomial if outlier = true ------------------------
  if (outlier & (out_prop > tol)) {
    Rcpp::NumericMatrix::Column zzcol = logprobs( Rcpp::_, logprobs.ncol() - 1);
    zzcol = dbetabinom_mu_rho_cpp(ocounts, osize, out_mean, out_disp, true) +
                                    log(out_prop);
  }

  // Log sum exponential trick of rows
  Rcpp::NumericVector ldensevec = logsumexp(logprobs);

  return ldensevec;
}

//' Objective function for the offspring.
//'
//' @inheritParams obj_offspring_vec
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
double obj_offspring(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                     int ploidy, Rcpp::NumericVector prob_geno,
                     double bias_val = 1, double seq_error = 0,
                     double od_param = 0,
                     bool outlier = false, double out_prop = 0.01,
                     double out_mean = 0.5, double out_disp = 1.0 / 3.0) {
  Rcpp::NumericVector ldensevec = obj_offspring_vec(ocounts, osize, ploidy, prob_geno,
                                                    bias_val, seq_error, od_param,
                                                    outlier, out_prop, out_mean, out_disp);
  return Rcpp::sum(ldensevec);
}

//' Just a reparameterization of \code{\link{obj_offspring}}.
//'
//' @inheritParams obj_offspring_vec
//' @param s Same as \code{exp(bias_val)} in \code{\link{obj_offspring}}.
//' @param ell We have \code{seq_error = expit(ell)} from \code{\link{obj_offspring}}.
//' @param r Same as \code{log((1.0 - od_param) / od_param)} from \code{\link{obj_offspring}}.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double obj_offspring_reparam(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                             int ploidy, Rcpp::NumericVector prob_geno,
                             double s, double ell,
                             double r) {
  double tol = 2.0 * DBL_EPSILON;

  double eps = expit(ell); // sequencing error rate
  double tau = 1.0 / (std::exp(r) + 1.0); // overidspersion parameter
  double d = std::exp(s); // bias parameter

  if (tau > (1.0 - tol)) {
    tau = 1.0 - 2.0 * tol;
  }
  return obj_offspring(ocounts, osize, ploidy, prob_geno, d, eps, tau, false, 0.01, 0.5, 1.0 / 3.0);
}

//' Same thing as \code{\link{obj_offspring}}, but each sample's log-density has a weight.
//'
//' This is mostly used in the EM algorithm.
//'
//' @inheritParams obj_offspring_vec
//' @param weight_vec A vector of numerics between 0 and 1. They don't have to sum to 1.
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
double obj_offspring_weights(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                             Rcpp::NumericVector weight_vec,
                             int ploidy, Rcpp::NumericVector prob_geno,
                             double bias_val = 1, double seq_error = 0,
                             double od_param = 0,
                             bool outlier = false, double out_prop = 0.01,
                             double out_mean = 0.5, double out_disp = 1.0 / 3.0) {
  // Check input --------------------------------------------------------------
  if (weight_vec.size() != ocounts.size()) {
    Rcpp::stop("weight_vec and ocounts should have the same size.");
  }
  for (int i = 0; i < weight_vec.size(); i++){
    if ((weight_vec(i) < 0) || (weight_vec(i) > 1)) {
      Rcpp::stop("weight_vec should all be between 0 and 1 (inclusive).");
    }
  }

  if (outlier) {
    Rcpp::warning("outlier = true in obj_offspring_weights.\nThis doesn't make much sense since you can't be using it for the EM.\nAt least not correctly.");
  }

  // Calculate log of densities for each observations -------------------------
  Rcpp::NumericVector ldensevec = obj_offspring_vec(ocounts, osize, ploidy, prob_geno,
                                                    bias_val, seq_error, od_param,
                                                    outlier, out_prop, out_mean, out_disp);
  // return weighted sum of log-densities --------------------------------------
  return Rcpp::sum(ldensevec * weight_vec);
}

//' Reparameterization of \code{\link{obj_offspring_weights}}.
//'
//' It doesn't make any sense to have outlier = true since the weights are just for the EM already.
//'
//' @inheritParams obj_offspring_weights
//' @inheritParams obj_offspring_reparam
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double obj_offspring_weights_reparam(Rcpp::NumericVector ocounts,
                                     Rcpp::NumericVector osize,
                                     Rcpp::NumericVector weight_vec,
                                     int ploidy, Rcpp::NumericVector prob_geno,
                                     double s, double ell,
                                     double r) {
  double eps = expit(ell); // sequencing error rate
  double tau = 1.0 / (std::exp(r) + 1.0); // over-dispersion parameter
  double d = std::exp(s); // bias parameter
  return obj_offspring_weights(ocounts, osize, weight_vec,
                               ploidy, prob_geno, d, eps,
                               tau, false, 0.01, 0.5, 1.0 / 3.0);
}
