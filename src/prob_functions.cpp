#include <Rcpp.h>

//' Returns the probability of seeing the reference allele after including
//' the mapping-bias and the sequencing-error.
//'
//' @param prob A numeric vector. Each element is the proportion of genomes that contain
//'     the reference allele. This should take on
//'     one of the values 0/K, 1/K, ... , K/K, where K is the ploidy of the individual.
//' @param bias The bias parameter. Should be greater than or equal to zero,
//'     though is typically less than 1.
//' @param seq_error The sequencing error rate. This should be between 0 and 1.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector pbias(Rcpp::NumericVector prob, double bias, double seq_error) {
  // check input ------------------------------------------

  for (int i = 0; i < prob.size(); i++){
    if (prob[i] < 0 || prob[i] > 1) {
      Rcpp::Rcout << "prob[" << i << "] = " << prob[i] << std::endl;
      Rcpp::stop("prob must be between 0 and 1.");
    }
  }
  if (bias < 0) {
    Rcpp::stop("bias must be greater than 0.");
  }
  if (seq_error < 0 || seq_error > 1) {
    Rcpp::stop("prob must be between 0 and 1.");
  }

  // Calculate the biased probability ---------------------
  Rcpp::NumericVector qvec = prob * (1 - seq_error) + (1 - prob) * seq_error;
  Rcpp::NumericVector xivec = qvec / (bias * (1 - qvec) + qvec);
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
Rcpp::NumericVector dbetabinom_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, double alpha_shape,
                             double beta_shape, bool return_log = false) {
  double tol = 2 * DBL_EPSILON; // tolerance from 0.

  // Check Input ------------------------------------------------
  if (alpha_shape < 0) {
    Rcpp::stop("alpha_shape must be greater than or equal to 0.");
  }
  if (beta_shape < 0) {
    Rcpp::stop("beta_shape must be greater than or equal to 0");
  }
  if (std::abs(beta_shape) <= tol && std::abs(alpha_shape) <= tol) {
    Rcpp::stop("beta_shape and alpha_shape cannot both be 0.");
  }
  if (x.size() != size.size()) {
    Rcpp::stop("x and size must have the same length.");
  }

  // Get ldense -------------------------------------------------
  Rcpp::NumericVector ldense(x.size());

  // RCpp::lbeta is for Rcpp::NumericVectors, R::lbeta is for doubles.
  if (std::abs(alpha_shape) > tol & std::abs(beta_shape) > tol) {
    ldense = Rcpp::lchoose(size, x) + Rcpp::lbeta(x + alpha_shape, size - x + beta_shape) -
      R::lbeta(alpha_shape, beta_shape);
  } else if (std::abs(alpha_shape) <= tol) {
    for (int i; i < ldense.size(); i++) {
      if (std::abs(x(i)) <= tol) {
        ldense(i) = 0;
      } else {
        ldense(i) = R_NegInf;
      }
    }
  } else if (std::abs(beta_shape) <= tol) {
    for (int i; i < ldense.size(); i++) {
      if (std::abs(x(i) - size(i)) <= tol) {
        ldense(i) = 0;
      } else {
        ldense(i) = R_NegInf;
      }
    }
  }

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
Rcpp::NumericVector dbetabinom_mu_rho_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, double mu,
                                    double rho, bool return_log = false) {
  double tol = 2 * DBL_EPSILON; // tolerance from 0.

  // Check input -------------------------------------------------------------
  if (mu < 0 || mu > 1) {
    Rcpp::stop("mu must be between 0 and 1 (inclusive).");
  }
  if (rho < tol || rho > 1 - tol) {
    Rcpp::stop("rho must be between 0 and 1 (not inclusive).");
  }

  // convert to shape parameters and call dbetabinom_cpp ---------------------
  double alpha_shape = mu * (1.0 - rho) / rho;
  double beta_shape  = (1.0 - mu) * (1.0 - rho) / rho;
  return dbetabinom_cpp(x, size, alpha_shape, beta_shape, return_log);
}
