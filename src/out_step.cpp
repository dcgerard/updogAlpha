// outlier objective functions and gradients

#include "updog.h"

//' Objective function of outlier part in EM step
//'
//' @param ocounts A vector of reads counts of reference allele.
//' @param osize A vector of total reads counts.
//' @param weight_vec The current probability of a point being an outlier.
//' @param out_mean The current mean of the beta.
//' @param out_disp The current overdispersion parameter of the beta.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double outlier_obj(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                   Rcpp::NumericVector weight_vec, double out_mean, double out_disp) {
  Rcpp::NumericVector ldensevec = dbetabinom_mu_rho_cpp(ocounts, osize, out_mean, out_disp, true);
  double obj = Rcpp::sum(weight_vec * ldensevec);
  return obj;
}


//' Derivative of overdispersion parameter with mean already calculated
//'
//' @inheritParams dbeta_dprop
//'
//' @author David Gerard
//'
//' @seealso dbeta_dtau
//'
// [[Rcpp::export]]
double dbeta_dtau_withxi(double x, double n, double xi, double tau){
  double h = (1.0 - tau) / tau;
  double dbdh = dbeta_dh(x, n, xi, h);
  double dhdtau = dh_dtau(tau);
  double deriv = dbdh * dhdtau;
  return deriv;
}

//' The gradient of \code{\link{outlier_obj}}.
//'
//' @inheritParams outlier_obj
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector outlier_grad(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                    Rcpp::NumericVector weight_vec, double out_mean, double out_disp) {

  Rcpp::NumericVector denom_vec = dbetabinom_mu_rho_cpp(ocounts, osize, out_mean, out_disp, false);

  Rcpp::NumericVector grad_vec(2);
  for (int i = 0; i < ocounts.size(); i++) {
    grad_vec(0) = grad_vec(0) + dbeta_dprop(ocounts(i), osize(i), out_mean, out_disp) * weight_vec(i) / denom_vec(i);
    grad_vec(1) = grad_vec(1) + dbeta_dtau_withxi(ocounts(i), osize(i), out_mean, out_disp) * weight_vec(i) / denom_vec(i);
  }

  return grad_vec;
}


