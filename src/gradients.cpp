////////
// This file contains functions to calculate the gradient
////////

#include "updog.h"
#include <math.h>

//' Derivative of mean of beta density
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
//' h = (1 - tau) / tau, so tau = 1 / (h + 1).
//'
//' @inheritParams dbeta_dprop
//' @param h A parameterization of the overdispersion parameter.
//'     This is just constrained to be greater than 0.
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


