#include "marlib_Rcpp.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNN(double alpha, NumericMatrix A, NumericMatrix B,
                      double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTN(double alpha, NumericMatrix A, NumericMatrix B,
                        double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::NOTRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmNT(double alpha, NumericMatrix A, NumericMatrix B,
                        double beta, NumericMatrix C) {
  marlib::dgemm(marlib::NOTRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::DenseMatrixT());
  return C;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_dgemmTT(double alpha, NumericMatrix A, NumericMatrix B,
                        double beta, NumericMatrix C) {
  marlib::dgemm(marlib::TRANS(), marlib::TRANS(), alpha, A, B, beta, C, marlib::DenseMatrixT(), marlib::DenseMatrixT());
  return C;
}
