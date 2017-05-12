// Functions for EM algorithm


#include "updog.h"

//' E step in EM algorithm for a single individual.
//'
//' @inheritParams obj_offspring
//' @inheritParams dbeta_deps
//'
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
Rcpp::NumericVector get_out_prop(Rcpp::NumericVector ocounts, Rcpp::NumericVector osize,
                                 int ploidy, int p1geno, int p2geno,
                                 double d, double eps, double tau, double out_prop,
                                 double out_mean, double out_disp) {
  // get density from good mixture -------------------------------------------------------
  Rcpp::NumericVector ldensevec = obj_offspring_vec(ocounts, osize, ploidy, p1geno, p2geno,
                                                    d, eps, tau, false, out_prop,
                                                    out_mean, out_disp) +
                                                      log(1.0 - out_mean);


  // get density from bad mixture --------------------------------------------------------
  Rcpp::NumericVector loutdensevec = dbetabinom_mu_rho_cpp(ocounts, osize, out_mean, out_disp, true) +
    log(out_mean);

  // return vector
  Rcpp::NumericVector max_el = Rcpp::pmax(ldensevec, loutdensevec);
  Rcpp::NumericVector denomvec = Rcpp::log(Rcpp::exp(ldensevec - max_el) + Rcpp::exp(loutdensevec - max_el)) + max_el;
  Rcpp::NumericVector theta_vec = Rcpp::exp(loutdensevec - denomvec);

  //Rcpp::NumericVector temp = Rcpp::exp(ldensevec - denomvec);
  //Rcpp::Rcout << theta_vec << std::endl
  //            << temp  << std::endl;
  return theta_vec;
}
