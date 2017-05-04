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

  // Get ldense -------------------------------------------------
  Rcpp::NumericVector ldense(x.size());

  // RCpp::lbeta is for Rcpp::NumericVectors, R::lbeta is for doubles.
  if ((std::abs(alpha_shape) > tol) & (std::abs(beta_shape) > tol)) {
    ldense = Rcpp::lchoose(size, x) + Rcpp::lbeta(x + alpha_shape, size - x + beta_shape) -
      R::lbeta(alpha_shape, beta_shape);
  } else if (std::abs(alpha_shape) <= tol) {
    for (int i = 0; i < ldense.size(); i++) {
      if (std::abs(x(i)) <= tol) {
        ldense(i) = 0;
      } else {
        ldense(i) = R_NegInf;
      }
    }
  } else if (std::abs(beta_shape) <= tol) {
    for (int i = 0; i < ldense.size(); i++) {
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
Rcpp::NumericVector dbetabinom_mu_rho_cpp(Rcpp::NumericVector x,
                                          Rcpp::NumericVector size, double mu,
                                          double rho, bool return_log = false) {
  double tol = 2 * DBL_EPSILON; // tolerance from 0.

  // Check input -------------------------------------------------------------
  if ((mu < 0) || (mu > 1)) {
    Rcpp::stop("mu must be between 0 and 1 (inclusive).");
  }
  if ((rho < -tol) || (rho > 1 - tol)) {
    Rcpp::stop("rho must be in [0, 1).");
  }

  // return binomial density if rho = 0 -----------------------
  if (rho < tol) {
    Rcpp::NumericVector dvec(x.size());
    for (int index = 0; index < x.size(); index++) {
      dvec(index) = R::dbinom(x(index), size(index),
           mu, (int)return_log);
    }
    return dvec;
  }

  // convert to shape parameters and call dbetabinom_cpp ---------------------
  double alpha_shape = mu * (1.0 - rho) / rho;
  double beta_shape  = (1.0 - mu) * (1.0 - rho) / rho;
  return dbetabinom_cpp(x, size, alpha_shape, beta_shape, return_log);
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



