// Small utility functions.

#include "updog.h"

//' Log-sum-exponential trick that I use all the time.
//'
//' @param xx A matrix whose rows are to be log-sum-exponentiated.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector logsumexp(Rcpp::NumericMatrix xx) {
  Rcpp::NumericVector lse_vec(xx.nrow());
  double maxrow;
  for (int r = 0; r < xx.nrow(); r++) {
    Rcpp::NumericMatrix::Row zzrow = xx(r, Rcpp::_); // reference rth row.
    maxrow = Rcpp::max(zzrow);
    lse_vec(r) = log(Rcpp::sum(Rcpp::exp(zzrow - maxrow))) + maxrow;
  }
  return lse_vec;
}
