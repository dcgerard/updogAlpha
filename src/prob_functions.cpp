////////
// Low level density and probability functions.
////////

#include "updog.h"

//' Returns the probability of seeing the reference allele after including
//' the mapping-bias and the sequencing-error.
//'
//' @param prob A numeric vector. Each element is the proportion of genomes that contain
//'     the reference allele. This should take on
//'     one of the values 0/K, 1/K, ... , K/K, where K is the ploidy of the individual.
//' @param bias The bias parameter. Should be greater than or equal to zero,
//'     though is typically less than 1. A 1 indicates no bias. A value less than one
//'     indicates a bias towards the reference allele. A value greater than 1 indicates
//'     a bias towards the non-reference allele.
//' @param seq_error The sequencing error rate. This should be between 0 and 1.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector pbias(Rcpp::NumericVector prob, double bias, double seq_error) {
  // Iterate through pbias_double
  Rcpp::NumericVector xivec(prob.size());
  for (int i = 0; i < prob.size(); i++){
      xivec(i) = pbias_double(prob(i), bias, seq_error);
  }
  return xivec;
}

//' A double version of \code{\link{pbias}}.
//'
//' @inheritParams pbias
//'
//' @seealso \code{\link{pbias}}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double pbias_double(double prob, double bias, double seq_error) {
  // check input -------------------------------------------
  if ((prob < 0) || (prob > 1)) {
    Rcpp::stop("prob must be between 0 and 1.");
  }
  if (bias < 0) {
    Rcpp::stop("bias must be greater than 0.");
  }
  if ((seq_error < 0) || (seq_error > 1)) {
    Rcpp::stop("prob must be between 0 and 1.");
  }

  double qvec = prob * (1.0 - seq_error) + (1.0 - prob) * seq_error;
  double xivec = qvec / (bias * (1.0 - qvec) + qvec);
  return xivec;
}

//' An Rcpp version of the \code{\link{dbetabinom}}.
//'
//' @param x A numeric vector of the observed counts.
//' @param size A numeric vector of the number of trials.
//' @param alpha_shape A numeric scalar of the alpha parameter underlying beta.
//' @param beta_shape A numeric scalar of the beta parameter of the underlying beta.
//' @param return_log A logical scalar. Should we return the log of the density (\code{TRUE})
//'     or not (\code{FALSE})?
//'
//' @export
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector dbetabinom_cpp(Rcpp::NumericVector x,
                                   Rcpp::NumericVector size,
                                   double alpha_shape,
                                   double beta_shape,
                                   bool return_log = false) {
  double tol = 2 * DBL_EPSILON; // tolerance from 0.

  // Check Input ------------------------------------------------
  if (alpha_shape < 0) {
    Rcpp::stop("alpha_shape must be greater than or equal to 0.");
  }
  if (beta_shape < 0) {
    Rcpp::stop("beta_shape must be greater than or equal to 0");
  }
  if ((std::abs(beta_shape) <= tol) && (std::abs(alpha_shape) <= tol)) {
    Rcpp::stop("beta_shape and alpha_shape cannot both be 0.");
  }
  if (x.size() != size.size()) {
    Rcpp::stop("x and size must have the same length.");
  }

  // call dbetabinom_mu_rho_cpp
  double mu = alpha_shape / (alpha_shape + beta_shape);
  double rho = 1.0 / (alpha_shape + beta_shape + 1.0);

  Rcpp::NumericVector ldense = dbetabinom_mu_rho_cpp(x, size, mu, rho, true);

  if (return_log) {
    return ldense;
  } else {
    return Rcpp::exp(ldense);
  }
}

//' An Rcpp version of \code{\link{dbetabinom_mu_rho}}.
//'
//' @inheritParams dbetabinom_cpp
//' @param mu The mean of the underlying beta.
//' @param rho The overdispersion parameter of the underlying beta.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector dbetabinom_mu_rho_cpp(Rcpp::NumericVector x,
                                          Rcpp::NumericVector size, double mu,
                                          double rho, bool return_log = false) {
  double tol = 2 * DBL_EPSILON; // tolerance from 0.

  // Check input -------------------------------------------------------------
  if (x.size() != size.size()) {
    Rcpp::stop("x and size must have the same length.");
  }

  // iterate through dbetabinom_mu_rho_cpp_double
  Rcpp::NumericVector ldense_vec(x.size());
  for (int i = 0; i < x.size(); i++) {
    ldense_vec(i) = dbetabinom_mu_rho_cpp_double(x(i), size(i), mu, rho, true);
  }

  // convert to shape parameters and call dbetabinom_cpp ---------------------
  if (return_log) {
    return ldense_vec;
  } else {
    return Rcpp::exp(ldense_vec);
  }
}

//' double version of betabinomial density.
//'
//' @param x The number of successes
//' @param size The number of draws
//' @param mu The mean proportion.
//' @param rho The overdispersion parameter.
//' @param return_log A logical. Should we return the log density
//'     (\code{TRUE}) or not (\code{FALSE})?
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double dbetabinom_mu_rho_cpp_double(double x, double size, double mu,
                                    double rho, bool return_log) {
  double tol = 2 * DBL_EPSILON; // tolerance from 0.

  // Check input -------------------------------------------------------------
  if (rho > (1.0 - tol)) {
    rho = 1.0 - 2 * tol;
  }
  if ((mu < 0) || (mu > 1)) {
    Rcpp::stop("mu must be between 0 and 1 (inclusive).");
  }
  if ((rho < 0) || (rho > 1 - tol)) {
    Rcpp::stop("rho must be in [0, 1).");
  }
  if ((x < 0) || (x > size)) {
    Rcpp::stop("x must be between 0 and size");
  }

  // return binomial density if rho is zero
  if (rho < tol) {
    return R::dbinom(x, size, mu, (int)return_log);
  }

  // return indicator function if mu is zero or 1
  if (mu < tol) {
    if ((x < tol) && return_log) {
      return 0.0;
    } else if ((x < tol) && !return_log) {
      return 1.0;
    } else if ((x > tol) && return_log) {
      return R_NegInf;
    } else {
      return 0.0;
    }
  }
  if (mu > (1 - tol)) {
    if ((std::abs(x - size) < tol) && return_log) {
      return 0.0;
    } else if ((std::abs(x - size) < tol) && !return_log) {
      return 1.0;
    } else if ((std::abs(x - size) > tol) && return_log) {
      return R_NegInf;
    } else {
      return 0.0;
    }
  }

  double alpha_shape = mu * (1.0 - rho) / rho;
  double beta_shape  = (1.0 - mu) * (1.0 - rho) / rho;

  double ldense = R::lchoose(size, x) +
    R::lgammafn(x + alpha_shape) +
    R::lgammafn(size - x + beta_shape) +
    R::lgammafn(alpha_shape + beta_shape) -
    R::lgammafn(size + alpha_shape + beta_shape) -
    R::lgammafn(alpha_shape) -
    R::lgammafn(beta_shape);

  if (return_log) {
    return ldense;
  } else {
    return exp(ldense);
  }
}

//' Vectorized version of \code{\link[stats]{dhyper}} for C++ implementation.
//'
//' @param x Vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
//' @param m The number of white balls in the urn.
//' @param n The number of black balls in the urn.
//' @param k The number of balls drawn from the urn.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector dhyper_cpp(Rcpp::NumericVector x,
                               int m, int n, int k) {
  Rcpp::NumericVector dense_vec(x.size());
  for (int i = 0; i < x.size(); i++) {
    dense_vec(i) = R::dhyper(x(i), m, n, k, 0); // zero to return density and not log
  }
  return dense_vec;
}



//' Rcpp function to get genotype probabilities assuming either
//' F1 population given parental genotypes, the allele frequencies assuming
//' Hardy-Weinberg equilibrium, or just a uniform distribution.
//'
//' @param ploidy The ploidy of the species.
//' @param model Do we assume the genotypes are distributed from an F1 population (\code{"f1"}),
//'     according to Hardy-Weinberg (\code{hw}), or uniformly (\code{"uniform"})?
//' @param p1geno The first parental genotype if \code{model = "f1"}.
//' @param p2geno The second parental genotype if \code{model = "f1"}.
//' @param allele_freq The allele-frequency if \code{model = "hw"}.
//'
//' @author David Gerard
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector get_prob_geno(int ploidy, std::string model, int p1geno, int p2geno, double allele_freq) {

  Rcpp::NumericVector prob_vec(ploidy + 1);

  if (model == "f1") {
    arma::Cube<double> q_array = get_q_array_cpp(ploidy);
    for (int i = 0; i <= ploidy; i++) {
      prob_vec(i) = q_array(p1geno, p2geno, i);
    }
  } else if (model == "hw") {
    for (int i = 0; i <= ploidy; i++) {
      prob_vec(i) = R::dbinom(i, ploidy, allele_freq, false);
    }
  } else if (model == "uniform") {
    for (int i = 0; i <= ploidy; i++) {
      prob_vec(i) = 1.0 / ((double)ploidy + 1.0);
    }
  } else {
    Rcpp::stop("model needs to be one of 'f1', 'hw', or 'uniform'");
  }
  return prob_vec;
}



//' Rcpp implementation of \code{\link{get_q_array}}.
//'
//' This function will return the segregation proabilities
//'
//' @inheritParams get_q_array
//'
//' @seealso \code{\link{get_q_array}}.
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
arma::Cube<double> get_q_array_cpp(int ploidy) {
  // check even and at least 2 --------------
  if (ploidy % 2 != 0) {
    Rcpp::stop("ploidy must be even.");
  }
  if (ploidy < 2) {
    Rcpp::stop("ploidy must be at least 2.");
  }

  arma::Cube<double> qarray(ploidy + 1, ploidy + 1, ploidy + 1);


  for (int oindex = 0; oindex < ploidy + 1; oindex++) {
    for (int p1index = 0; p1index < ploidy + 1; p1index++) {
      for (int p2index = 0; p2index < ploidy + 1; p2index++) {
        if (p1index + p2index < oindex) {
          qarray(p1index, p2index, oindex) = 0.0;
        } else {
          int minval = std::max(0, oindex - p2index);
          int maxval = std::min(ploidy / 2, p1index);

          if (minval > maxval) {
            int trash_int = minval;
            minval = maxval;
            maxval = trash_int;
          }

          Rcpp::NumericVector aseq(maxval - minval + 1);
          for (int i = 0; i < aseq.size(); i++) {
            aseq(i) = minval + i;
          }
          Rcpp::NumericVector p1prob = dhyper_cpp(aseq,
                                                  p1index,
                                                  ploidy - p1index,
                                                  ploidy / 2);
          Rcpp::NumericVector p2prob = dhyper_cpp(oindex - aseq,
                                                  p2index,
                                                  ploidy - p2index,
                                                  ploidy / 2);
          qarray(p1index, p2index, oindex) = Rcpp::sum(p1prob * p2prob);
        }
      }
    }
  }
  return qarray;
}
