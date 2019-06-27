#include "marlib_Rcpp.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int C_poi_rightbound(double lambda, double eps) {
  return marlib::poi::rightbound(lambda, eps);
}

//' @export
// [[Rcpp::export]]
List C_poi_pmf(double lambda, int left, int right) {
  NumericVector prob(right-left+1);
  double weight = marlib::poi::pmf(lambda, left, right, prob);
  return List::create(Named("prob")=prob, Named("weight")=weight);
}

//' @export
// [[Rcpp::export]]
List C_poi_cpmf(double lambda, int left, int right) {
  NumericVector prob(right-left+1);
  NumericVector cprob(right-left+1);
  double weight = marlib::poi::cpmf(lambda, left, right, prob, cprob);
  return List::create(Named("prob")=prob, Named("cprob")=cprob, Named("weight")=weight);
}
