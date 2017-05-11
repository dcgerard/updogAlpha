////////
// This file contains functions to calculate the gradient
////////

#include "updog.h"
#include <math.h>

//' Derivative of mean of beta density.
//'
//' This calculates d/dxi BB(x|n, xi, tau).
//'
//' @param x The observed counts of reference allele.
//' @param n The total number of observed counts.
//' @param xi The mean proportion.
//' @param tau The overdispersion parameter.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dbeta_dprop(double x, double n, double xi, double tau) {
  // check input -----------------------------------------------
  if ((x > n) || (x < 0)) {
    Rcpp::stop("x must be between 0 and n");
  }
  if (n < 1) {
    Rcpp::stop("n must be greater than or equal to 1.");
  }
  if ((xi < 0) || (xi > 1)) {
    Rcpp::stop("xi must be between 0 and 1");
  }
  if ((tau < 0) || (tau > 1)) {
    Rcpp::stop("tau must be between 0 and 1");
  }

  double ldense = dbetabinom_mu_rho_cpp_double(x, n, xi, tau, true);

  // Compute all needed gamma and digamma functions ------------
  double h     = (1 - tau) / tau;
  double comp1 = h * R::digamma(x + xi * h); // in numer
  double comp2 = -1.0 * h * R::digamma(n - x + (1 - xi) * h); // in numer
  double comp3 = h * R::digamma(xi * h); // in denom
  double comp4 = -1.0 * h * R::digamma((1.0 - xi) * h); // in denom
  double deriv = (comp1 + comp2 - comp3 - comp4) * std::exp(ldense);

  return deriv;
}

//' Derivative of overdispersion parameter of beta density.
//'
//' This calculates d/dh of BB(x|n, xi, h), where h = (1 - tau) / tau, so tau = 1 / (h + 1).
//'
//' @inheritParams dbeta_dprop
//' @param h A parameterization of the overdispersion parameter.
//'     This is just constrained to be greater than 0.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dbeta_dh(double x, double n, double xi, double h) {

  // Check input
  if (h < 0) {
    Rcpp::stop("h must be greater than 0");
  }

  // get density ------------------------------------------
  double rho = 1.0 / (h + 1.0);
  double ldense = dbetabinom_mu_rho_cpp_double(x, n, xi, rho, true);

  // get six components -----------------------------------
  double comp1 =  xi * R::digamma(x + xi * h);
  double comp2 = (1.0 - xi) * R::digamma(n - x + (1.0 - xi) * h);
  double comp3 = R::digamma(h);
  double comp4 = R::digamma(n + h);
  double comp5 = xi * R::digamma(xi * h);
  double comp6 = (1.0 - xi) * R::digamma((1.0 - xi) * h);
  double deriv = (comp1 + comp2 + comp3 - comp4 - comp5 - comp6) * std::exp(ldense);

  return deriv;
}


//' Returns derivative to f / (d * (1 - f) + f)
//'
//' @param f The sequencing-error adjusted probability of the refernence allele.
//' @param d The bias term. Must be greater than 0.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dxi_df(double d, double f) {
  // Check input -----------------------------------------
  if ((f < 0) || (f > 1)) {
    Rcpp::stop("f must be between 0 and 1");
  }
  if (d < 0) {
    Rcpp::stop("d must be greater than 0");
  }

  double deriv = d / std::pow(d * (1.0 - f) + f, 2.0);
  return deriv;
}

//' Returns derivative of p(1 - eps) + (1 - p) * eps
//'
//' @param p The proportion of genomes whose allele is A.
//' @param eps The sequencing error rate.
//'
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double df_deps(double eps, double p) {
  if ((eps < 0) || (eps > 1)) {
    Rcpp::stop("eps needs to be between 0 and 1 (inclusive).");
  }
  if ((p < 0) || (p > 1)) {
    Rcpp::stop("p needs to be between 0 and 1 (inclusive).");
  }
  return 1.0 - 2.0 * p;
}

//' Derivative of exp(ell) / (1 + exp(ell))
//'
//' exp(ell) / (1 + exp(ell)) is the sequencing error rate and ell is a parameterization to
//' that is unconstrained.
//'
//' @param ell A double.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double deps_dell(double ell) {
  double deriv = 1.0 / std::pow(std::exp(ell / -2.0) + std::exp(ell / 2.0), 2.0);
  return(deriv);
}



//' Derivative of beta density w.r.t. unconstrained parameterization of sequencing error.
//'
//' I use the chain rule here. \code{\link{dbeta_dprop}} * \code{\link{dxi_df}} *
//' \code{\link{df_deps}} * \code{\link{deps_dell}}.
//'
//' @param x The number of counts of reference allele.
//' @param n The number of counts of reads.
//' @param d The sequencing bias parameter.
//' @param ell Expit of eps. We have eps = exp(ell) / (1 + exp(ell)), or ell = log(eps / (1 - eps))
//' @param p The proportion of genome that is the reference allele.
//' @param h the overdisperion parameter. tau = 1 / (h + 1) or h = (1 - tau) / tau
//'
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dbeta_dl(double x, double n, double d, double ell, double p, double h) {

  // intermediate parameters --------------------------------------------------
  double tau = 1.0 / (h + 1.0);
  double eps = expit(ell);
  double xi = pbias_double(p, d, eps);
  double f = (1.0 - p) * eps + p * (1.0 - eps);

  // Derivatives --------------------------------------------------------------
  double dbdxi  = dbeta_dprop(x, n, xi, tau);
  double dxidf  = dxi_df(d, f);
  double dfdeps = df_deps(eps, p);
  double depsdl = deps_dell(ell);
  double deriv = dbdxi * dxidf * dfdeps * depsdl;
  return deriv;
}

//' Derivative w.r.t. d of xi(d, f) = f / (d * (1 - f) + f)
//'
//' @inheritParams dxi_df
//'
//' @author David Gerard
//'
//'
// [[Rcpp::export]]
double dxi_dd(double d, double f) {
  double deriv = -1.0 * f * (1.0 - f) / std::pow(d * (1.0 - f) + f, 2);
  return deriv;
}

//' Derivative of betabinomial density w.r.t. bias parameter.
//'
//' Uses chain rule with \code{\link{dbeta_dprop}} * \code{\link{dxi_dd}}.
//'
//' @inheritParams dbeta_dl
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dbeta_dd(double x, double n, double d, double ell, double p, double h) {

  // intermediate parameters --------------------------------------------------
  double tau = 1.0 / (h + 1.0);
  double eps = expit(ell); // sequencing error
  double xi = pbias_double(p, d, eps); // adjusted prob of A
  double f = (1.0 - p) * eps + p * (1.0 - eps);

  // derivatives --------------------------------------------------------------
  double dbdxi = dbeta_dprop(x, n, xi, tau);
  double dxidd = dxi_dd(d, f);
  double deriv = dbdxi * dxidd;
  return deriv;
}

//' Same as \code{\link{dbeta_dh}}, but with same inputs as \code{\link{dbeta_dd}}
//' and \code{\link{dbeta_dl}}.
//'
//' @inheritParams dbeta_dl
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dbeta_dh_ell(double x, double n, double d, double ell, double p, double h) {
  // intermediate parameters ---------------------------------------------------
  double eps = expit(ell);
  double xi = pbias_double(p, d, eps);

  // derivative ----------------------------------------------------------------
  double deriv = dbeta_dh(x, n, xi, h);
  return deriv;
}


//' Gradient of \code{\link{obj_offspring}} for each individual.
//'
//' @inheritParams obj_offspring
//' @param d The bias term
//' @param ell The logit of the sequencing error rate
//' @param h We have h = (1 - tau) / tau
//'
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix grad_offspring_mat(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                                       int ploidy, int p1geno, int p2geno,
                                       double d, double ell,
                                       double h) {
  double tol = 2.0 * DBL_EPSILON;

  // Get ancillary parameters ------------------------------------------------------
  double tau = 1.0 / (h + 1.0); // the overdispersion parameter between 0 and 1
  double eps = expit(ell); // the sequencing error rate

  // Get the log of the denominator for each individual ----------------------------
  Rcpp::NumericVector ldenom_vec = obj_offspring_vec(ocounts, osize,
                                                     ploidy, p1geno, p2geno,
                                                     d, eps, tau, false, 0, 1.0 / 2.0, 1.0 / 3.0);

  // Get possible probabilities ----------------------------------------------------
  Rcpp::NumericVector probs(ploidy + 1);
  for (int i = 0; i < ploidy + 1; i++) {
    probs(i) = (double)i / ploidy;
  }

  // Get segregation probabilities --- probably faster to have this as argument ----
  arma::Cube<double> qarray = get_q_array_cpp(ploidy);

  Rcpp::NumericMatrix grad_ind(ocounts.size(), 3); // goes d, ell, h
  for (int i = 0; i < ocounts.size(); i++) {
    for (int j = 0; j < ploidy + 1; j++) {
      if (qarray(p1geno, p2geno, j) > tol) {
        grad_ind(i, 0) = grad_ind(i, 0) + dbeta_dd(ocounts(i), osize(i), d, ell, probs(j), h) *
          qarray(p1geno, p2geno, j);
        grad_ind(i, 1) = grad_ind(i, 1) + dbeta_dl(ocounts(i), osize(i), d, ell, probs(j), h) *
          qarray(p1geno, p2geno, j);
        grad_ind(i, 2) = grad_ind(i, 2) + dbeta_dh_ell(ocounts(i), osize(i), d, ell, probs(j), h) *
          qarray(p1geno, p2geno, j);
      }
    }
    grad_ind(i, 0) = grad_ind(i, 0) * std::exp(-1.0 * ldenom_vec(i));
    grad_ind(i, 1) = grad_ind(i, 1) * std::exp(-1.0 * ldenom_vec(i));
    grad_ind(i, 2) = grad_ind(i, 2) * std::exp(-1.0 * ldenom_vec(i));
  }

  return grad_ind;
}

//' Gradient of \code{\link{obj_offspring_reparam}}.
//'
//' @inheritParams obj_offspring
//' @inheritParams grad_offspring_mat
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector grad_offspring(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                                   int ploidy, int p1geno, int p2geno,
                                   double d, double ell,
                                   double h) {
  Rcpp::NumericMatrix gb_mat = grad_offspring_mat(ocounts, osize, ploidy, p1geno, p2geno,
                                                  d, ell, h);
  Rcpp::NumericVector grad = colSums_cpp(gb_mat);
  return(grad);
}

//' Gradient of \code{\link{obj_offspring_weights_reparam}}
//'
//' @inheritParams grad_offspring
//' @param weight_vec A vector of weights between 0 and 1 (do not need to add up to 1).
//'
// [[Rcpp::export]]
Rcpp::NumericVector grad_offspring_weights(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                                           Rcpp::NumericVector weight_vec,
                                           int ploidy, int p1geno, int p2geno,
                                           double d, double ell,
                                           double h) {
  // Check input --------------------------------------------------------------
  if (weight_vec.size() != ocounts.size()) {
    Rcpp::stop("weight_vec and ocounts should have the same size.");
  }
  for (int i = 0; i < weight_vec.size(); i++){
    if ((weight_vec(i) < 0) || (weight_vec(i) > 1)) {
      Rcpp::stop("weight_vec should all be between 0 and 1 (inclusive).");
    }
  }
  Rcpp::NumericMatrix gb_mat = grad_offspring_mat(ocounts, osize, ploidy, p1geno, p2geno,
                                                  d, ell, h);
  Rcpp::NumericVector grad(3);
  Rcpp::NumericMatrix::Column zzcol = gb_mat(Rcpp::_, 0);
  grad(0) = Rcpp::sum(zzcol * weight_vec);
  zzcol = gb_mat(Rcpp::_, 1);
  grad(1) = Rcpp::sum(zzcol * weight_vec);
  zzcol = gb_mat(Rcpp::_, 2);
  grad(2) = Rcpp::sum(zzcol * weight_vec);
  return(grad);
}


