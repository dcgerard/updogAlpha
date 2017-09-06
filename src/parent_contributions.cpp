// objective and gradient contributions from parental sequences.

#include "updog.h"

//' Parent contribution to objective function.
//'
//' For each parental read counts you have, add this to
//' \code{{\link{obj_offspring}}} to get the objective function.
//'
//' @inheritParams obj_offspring_vec
//' @param pcounts The number of counts of the reference allele we observe in
//'     a parent.
//' @param psize The total number of reads we observe in a parent.
//' @param pgeno The genotype of the parent.
//' @param weight A double that multiplies to the log-likelihood (objective).
//'     Should only be used when \code{outlier = FALSE} during the M
//'     step of the EM algorithm.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double obj_parent(double pcounts, double psize,
                  int ploidy, int pgeno,
                  double bias_val = 1, double seq_error = 0,
                  double od_param = 0,
                  bool outlier = false, double out_prop = 0.01,
                  double out_mean = 0.5, double out_disp = 1.0 / 3.0,
                  double weight = 1.0) {
  // Check input ---------------------------------------------------------
  double tol = 2 * DBL_EPSILON;
  if ((weight < -1.0 * tol) | (weight > 1.0 + tol)) {
    Rcpp::stop("weight must be between 0 and 1.");
  }
  if ((weight < 1.0 - tol) & outlier) {
    Rcpp::stop("having weight < 1.0 and outlier = true makes no sense.\nOnly use the weight parameter when doing the em.");
  }
  if ((bias_val) < 0.0) {
    Rcpp::stop("bias_val needs to be greater than or equal to 0.");
  }
  if ((seq_error) < 0.0 | (seq_error) > 1.0) {
    Rcpp::stop("seq_error needs to be between 0 and 1.");
  }
  if ((od_param) < 0.0 | (od_param > 1.0)) {
    Rcpp::stop("od_param needs to be between 0 and 1.");
  }
  if ((pcounts > psize) | (pcounts < 0)) {
    Rcpp::stop("pcounts needs to be between 0 and psize");
  }

  // Calculate beta-binomial mixture density -------------------------------


  double bbmean = pbias_double((double)pgeno / (double)ploidy, bias_val, seq_error);
  double good_obj = dbetabinom_mu_rho_cpp_double(pcounts, psize, bbmean, od_param, true);

  double pobj;
  if (outlier & (out_prop > tol)) {
    good_obj = good_obj + log(1 - out_prop);
    double out_obj = dbetabinom_mu_rho_cpp_double(pcounts, psize, out_mean, out_disp, true) +
      log(out_prop);
    double max_obj = std::max(good_obj, out_obj);
    pobj = max_obj + log(std::exp(good_obj - max_obj) + std::exp(out_obj - max_obj));
  } else {
    pobj = good_obj;
  }

  return weight * pobj;
}

//' A reparameterization of \code{\link{obj_parent}}
//'
//' @inheritParams obj_offspring_reparam
//' @inheritParams obj_offspring_vec
//' @inheritParams obj_parent
//' @param weight A double between 0 and 1. This should only not be 1 when
//'     doing the em algorithm and the weight is the probabiility that
//'     a point being ok.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double obj_parent_reparam(double pcounts, double psize, int ploidy,
                          int pgeno, double s, double ell, double r,
                          double weight = 1.0,
                          bool outlier = false, double out_prop = 0.01) {
  double tol = 2.0 * DBL_EPSILON;
  if ((weight < -1.0 * tol) | (weight > 1.0 + tol)) {
    Rcpp::stop("weight must be between 0 and 1.");
  }

  double eps = expit(ell); // sequencing error rate
  double tau = 1.0 / (std::exp(r) + 1.0); // overidspersion parameter
  double d = std::exp(s); // bias parameter

  if (tau > (1.0 - tol)) {
    tau = 1.0 - 2.0 * tol;
  }
  return weight * obj_parent(pcounts, psize, ploidy, pgeno, d, eps, tau, outlier, out_prop, 0.5, 1.0 / 3.0);
}

//////////////////////////////////////////////////
// Gradients
/////////////////////////////////////////////////

//' Gradient of \code{\link{obj_parent}} (when \code{outlier = FALSE}) with respect to
//' \code{bias_val}, \code{seq_error}, and \code{od_param}. Basically, this just calculates
//' derivatives of the log beta-binomial density.
//'
//' @inheritParams obj_parent
//' @inheritParams obj_offspring_vec
//'
//' @return A \code{NumericVector} of length three that contains the partial derivatives of
//'     \code{bias_val}, \code{seq_error}, and \code{od_param}, in that order.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector grad_parent(double pcounts, double psize,
                                int ploidy, int pgeno,
                                double bias_val = 1, double seq_error = 0,
                                double od_param = 0,
                                double weight = 1.0) {

  // Check Input ------------------------------------------------------------
  double tol = 2.0 * DBL_EPSILON;
  if ((bias_val) < 0.0) {
    Rcpp::stop("bias_val needs to be greater than or equal to 0.");
  }
  if ((seq_error) < 0.0 | (seq_error) > 1.0) {
    Rcpp::stop("seq_error needs to be between 0 and 1.");
  }
  if ((od_param) < 0.0 | (od_param > 1.0)) {
    Rcpp::stop("od_param needs to be between 0 and 1.");
  }
  if ((pcounts > psize) | (pcounts < 0)) {
    Rcpp::stop("pcounts needs to be between 0 and psize");
  }
  if ((weight < 0.0) | (weight > 1.0)) {
    Rcpp::stop("weight needs to be between 0 and 1.");
  }

  // Calculate gradient ---------------------------------------------------
  double p         = (double)pgeno / (double)ploidy; // naive prob
  double mu        = pbias_double(p, bias_val, seq_error); // mean of bb
  double fval      = dbetabinom_mu_rho_cpp_double(pcounts, psize, mu, od_param, false);
  double grad_bias = weight * dbeta_dd(pcounts, psize, bias_val, seq_error, p, od_param) / fval;
  double grad_seq  = weight * dbeta_deps(pcounts, psize, bias_val, seq_error, p, od_param) / fval;
  double grad_od   = weight * dbeta_dtau(pcounts, psize, bias_val, seq_error, p, od_param) / fval;

  Rcpp::NumericVector grad(3);
  grad(0) = grad_bias;
  grad(1) = grad_seq;
  grad(2) = grad_od;

  return(grad);
}


//' The gradient of \code{\link{obj_parent_reparam}} with respect to \code{s}, \code{ell},
//' and \code{r}, which recall are reparameterizations of the bias, seqeuncing error, and
//' overdispersion, respectively.
//'
//' @inheritParams obj_offspring_reparam
//' @inheritParams obj_offspring_vec
//' @inheritParams obj_parent
//'
//' @return A \code{NumericVector} of length three with the partial derivatives of \code{s},
//'     \code{ell}, and \code{r} in that order.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector grad_parent_reparam(double pcounts, double psize, int ploidy,
                                        int pgeno, double s, double ell, double r,
                                        double weight = 1.0) {

  double p        = (double)pgeno / (double)ploidy; // original prob.
  double h        = std::exp(r); // another parameterization of the overdispersion.
  double d        = std::exp(s); // another parameterization of the bias.
  double tau      = 1.0 / (h + 1.0); // original od parameterization.
  double eps      = expit(ell); // sequencing error rate.
  double mu       = pbias_double(p, d, eps); // mean of bb.
  double fval     = dbetabinom_mu_rho_cpp_double(pcounts, psize, mu, tau, false);
  double grad_s   = weight * dbeta_ds(pcounts, psize, s, ell, p, h) / fval;
  double grad_ell = weight * dbeta_dl(pcounts, psize, d, ell, p, h) / fval;
  double grad_r   = weight* dbeta_dr_ell(pcounts, psize, d, ell, p, r) / fval;

  Rcpp::NumericVector grad(3);
  grad(0) = grad_s;
  grad(1) = grad_ell;
  grad(2) = grad_r;

  return(grad);
}

//' Gradient function for \code{\link{dbetabinom_mu_rho_cpp_double}}.
//'
//' @inheritParams obj_parent
//' @param out_mean The mean of the BB.
//' @param out_disp The overdispersion parameter of the BB.
//' @param weight The probability a point is an outlier.
//'
//' @return A vector of length two. The first element of which is the partial derivative
//'     of the outlier mean. The second element of which is the partial
//'     derivative of the outlier overdispersion parameter.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector out_grad_parent(double pcounts, double psize,
                                    double out_mean, double out_disp,
                                    double weight) {

  Rcpp::NumericVector grad(2);

  double fval = dbetabinom_mu_rho_cpp_double(pcounts, psize, out_mean, out_disp, false);
  grad(0) = weight * dbeta_dprop(pcounts, psize, out_mean, out_disp) / fval;
  double h = (1.0 - out_disp) / out_disp;
  grad(1) = weight * dbeta_dh(pcounts, psize, out_mean, h) * dh_dtau(out_disp) / fval;

  return grad;
}


//' E-step for the parents.
//'
//' @inheritParams obj_parent
//' @inheritParams obj_offspring_vec
//' @param d Same as \code{bias_val} in \code{\link{obj_parent}}.
//' @param eps Same as \code{seq_error} in \code{\link{obj_parent}}.
//' @param tau Same as \code{od_param} in \code{\link{obj_parent}}.
//'
//' @author David Gerard
//'
//'
// [[Rcpp::export]]
double get_parent_outprop(double pcounts, double psize,
                          int ploidy, int pgeno,
                          double d, double eps, double tau,
                          double out_prop, double out_mean,
                          double out_disp) {
  double tol  = 2.0 * DBL_EPSILON;
  double prob = (double)pgeno / (double)ploidy;
  double xi   = pbias_double(prob, d, eps);

  double theta; // The weight to return.
  if (out_prop < tol) {
    theta = 0.0;
  } else {
    double lgood   = dbetabinom_mu_rho_cpp_double(pcounts, psize, xi, tau, true) + std::log(1 - out_prop);
    double lout    = dbetabinom_mu_rho_cpp_double(pcounts, psize, out_mean, out_disp, true) + log(out_prop);
    double max_val = std::max(lgood, lout);
    double ldenom  = std::log(std::exp(lgood - max_val) + std::exp(lout - max_val)) + max_val;
    theta          = std::exp(lout - ldenom);
  }

  return theta;
}



