#include "marlib_Rcpp.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double C_ddot(NumericVector x, NumericVector y) {
  return marlib::ddot(x,y);
}

//' @export
// [[Rcpp::export]]
double C_ddot2(S4 x, NumericVector y) {
  return marlib::ddot(x,y);
}

//' @export
// [[Rcpp::export]]
double C_dasum(NumericVector x) {
  return marlib::dasum(x);
}

//' @export
// [[Rcpp::export]]
double C_idamax(NumericVector x) {
  return marlib::idamax(x);
}

//' @export
// [[Rcpp::export]]
NumericVector C_dcopy(NumericVector x, NumericVector y) {
  marlib::dcopy(x, y);
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dcopy2(S4 x, NumericVector y) {
  marlib::dcopy(x, y);
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dscal(double alpha, NumericVector x) {
  marlib::dscal(alpha, x);
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_daxpy(double alpha, NumericVector x, NumericVector y) {
  marlib::daxpy(alpha, x, y);
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dfill(NumericVector x, double v) {
  marlib::dfill(x, v, marlib::ArrayT());
  return x;
}
