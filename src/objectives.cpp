////////
// This file contains updog objective functions
////////

#include "updog.h"

//' Objective function for the offspring.
//'
//' @param ocounts A vector of integers. The observed counts of the refernce
//'     allele for each individual.
//' @param osize A vector of integers. The observed number of reads for each
//'     individuals.
//' @param ploidy An integer. The ploidy of the species. This is assumed
//'     to be the same for each individual.
//' @param p1geno The first parental genotype. The number of copies of the
//'     reference allele the first parent has.
//' @param p2geno The second parental genotype. The number of copies of the
//'     reference allele the second parent has.
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
// [[Rcpp::export]]
double obj_offspring(Rcpp::IntegerVector ocounts, Rcpp::IntegerVector osize,
                     int ploidy, int p1geno, int p2geno,
                     double bias_val = 1, double seq_error = 0,
                     double od_param = 0,
                     bool outlier = false, double out_prop = 0.01,
                     double out_mean = 0.5, double out_disp = 1.0 / 3.0) {
  double tol = 2 * DBL_EPSILON;

  // check input ---------------------------------------------------
  if (ocounts.size() != osize.size()) {
    Rcpp::stop("ocounts and osize must have the same length.");
  }
  if ((p1geno < 0) | (p1geno > ploidy)) {
    Rcpp::stop("p1geno must be between 0 and ploidy");
  }
  if ((p2geno < 0) + (p2geno > ploidy)) {
    Rcpp::stop("p2geno must be between 0 and ploidy");
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
  Rcpp::NumericVector prob(ploidy + 1);
  for (int i = 0; i < ploidy + 1; i++) {
    prob(i) = (double)i / ploidy;
  }
  prob = pbias(prob, bias_val, seq_error);

  // calculate segregation probabilities
  // might be faster to compute this just once and pass it.
  arma::Cube<double> qarray = get_q_array_cpp(ploidy);

  // Calculate probabilities for OK points --------------------------
  // only need last column is outlier = true.

  // count how many non-zero values in qarray(p1geno, p2geno, ) there are
  int colnum = 0;
  for (int i = 0; i < ploidy + 1; i++) {
    if (qarray(p1geno, p2geno, i) > tol) {
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
    if (qarray(p1geno, p2geno, i) > tol) {
      Rcpp::NumericMatrix::Column zzcol = logprobs(Rcpp::_, colindex); // reference to colindexth column
      zzcol = dbetabinom_mu_rho_cpp((Rcpp::NumericVector)ocounts,
                                    (Rcpp::NumericVector)osize,
                                    prob(i), od_param, true) +
                                      log(qarray(p1geno, p2geno, i));
      if (outlier) {
        zzcol = zzcol + log(1 - out_prop);
      }
      colindex++;
    }
  }

  // Another beta binomial if outlier = true ------------------------
  if (outlier & (out_prop > tol)) {
    Rcpp::NumericMatrix::Column zzcol = logprobs( Rcpp::_, logprobs.ncol() - 1);
    zzcol = dbetabinom_mu_rho_cpp((Rcpp::NumericVector)ocounts,
                                  (Rcpp::NumericVector)osize,
                                  out_mean, out_disp, true) +
                                    log(out_prop);
  }

  // Log sum exponential trick of rows
  double ldense = Rcpp::sum(logsumexp(logprobs));

  return ldense;
}


