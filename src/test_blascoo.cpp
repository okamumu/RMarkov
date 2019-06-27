#include "marlib_Rcpp.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvN_coo(double alpha, S4 A, NumericVector x,
                      double beta, NumericVector y) {
  marlib::dgemv(marlib::NOTRANS(), alpha, A, x, beta, y, marlib::COOMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
NumericVector C_dgemvT_coo(double alpha, S4 A, NumericVector x,
                           double beta, NumericVector y) {
  marlib::dgemv(marlib::TRANS(), alpha, A, x, beta, y, marlib::COOMatrixT());
  return y;
}

//' @export
// [[Rcpp::export]]
S4 C_dgerN_coo(double alpha, NumericVector x, NumericVector y, S4 A) {
  marlib::dger(marlib::NOTRANS(), alpha, x, y, A, marlib::COOMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
S4 C_dgerT_coo(double alpha, NumericVector x, NumericVector y, S4 A) {
  marlib::dger(marlib::TRANS(), alpha, x, y, A, marlib::COOMatrixT());
  return A;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN_coo(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::COOMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN_coo(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::COOMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT_coo(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::COOMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT_coo(double alpha, S4 A, NumericMatrix B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::COOMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN2_coo(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::COOMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN2_coo(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::COOMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT2_coo(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::COOMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT2_coo(double alpha, NumericMatrix A, S4 B, double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::COOMatrixT());
  return C;
}
