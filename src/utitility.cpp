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

//' Gets all possible binomial probabilities for a given ploidy, bias term, and sequencing
//' error rate.
//'
//' @param ploidy The ploidy of the species.
//' @param bias_val The bias parameter. Should be greater than 0. A value of 1 means no bias.
//' @param seq_error The sequencing error rate.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector get_pvec(int ploidy, double bias_val, double seq_error) {
  Rcpp::NumericVector prob(ploidy + 1);
  for (int i = 0; i < ploidy + 1; i++) {
    prob(i) = (double)i / ploidy;
  }
  prob = pbias(prob, bias_val, seq_error);
  return prob;
}


//' Stupid implementation of colSums because I guess not implemented in Rcpp sugar.
//'
//' @param x A NumericMatrix.
//'
//' @author David Gerard
//'
//'
// [[Rcpp::export]]
Rcpp::NumericVector colSums_cpp(Rcpp::NumericMatrix x) {
  Rcpp::NumericVector sumvec(x.ncol());
  for (int i = 0; i < x.nrow(); i++) {
    for (int j = 0; j < x.ncol(); j++) {
      sumvec(j) = sumvec(j) + x(i, j);
    }
  }
  return(sumvec);
}

//' The expit function.
//'
//' @param x A double.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double expit(double x) {
  return std::exp(x) / (1 + std::exp(x));
}
