#include <marlib.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector C_gsstep_vecN_dense(double alpha, NumericMatrix A,
                                  double sigma, double omega, NumericVector b, NumericVector x0) {
  NumericVector x = clone(x0);
  marlib::gsstep(marlib::NOTRANS(), alpha, A, sigma, omega, b, x,
                 marlib::DenseMatrixT(), marlib::ArrayT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_gsstep_vecT_dense(double alpha, NumericMatrix A,
                                  double sigma, double omega, NumericVector b, NumericVector x0) {
  NumericVector x = clone(x0);
  marlib::gsstep(marlib::TRANS(), alpha, A, sigma, omega, b, x,
                 marlib::DenseMatrixT(), marlib::ArrayT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_gsstep_vecN_csr(double alpha, S4 A,
                                  double sigma, double omega, NumericVector b, NumericVector x0) {
  NumericVector x = clone(x0);
  marlib::gsstep(marlib::NOTRANS(), alpha, A, sigma, omega, b, x,
                 marlib::CSRMatrixT(), marlib::ArrayT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericVector C_gsstep_vecT_csc(double alpha, S4 A,
                                double sigma, double omega, NumericVector b, NumericVector x0) {
  NumericVector x = clone(x0);
  marlib::gsstep(marlib::TRANS(), alpha, A, sigma, omega, b, x,
                 marlib::CSCMatrixT(), marlib::ArrayT());
  return x;
}

/// gs mat

//' @export
// [[Rcpp::export]]
NumericMatrix C_gsstep_matN_dense(double alpha, NumericMatrix A,
                                  double sigma, double omega, NumericMatrix b, NumericMatrix x0) {
  NumericMatrix x = clone(x0);
  marlib::gsstep(marlib::NOTRANS(), alpha, A, sigma, omega, b, x,
                 marlib::DenseMatrixT(), marlib::DenseMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_gsstep_matT_dense(double alpha, NumericMatrix A,
                                  double sigma, double omega, NumericMatrix b, NumericMatrix x0) {
  NumericMatrix x = clone(x0);
  marlib::gsstep(marlib::TRANS(), alpha, A, sigma, omega, b, x,
                 marlib::DenseMatrixT(), marlib::DenseMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_gsstep_matN_csr(double alpha, S4 A,
                                  double sigma, double omega, NumericMatrix b, NumericMatrix x0) {
  NumericMatrix x = clone(x0);
  marlib::gsstep(marlib::NOTRANS(), alpha, A, sigma, omega, b, x,
                 marlib::CSRMatrixT(), marlib::DenseMatrixT());
  return x;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_gsstep_matT_csc(double alpha, S4 A,
                                  double sigma, double omega, NumericMatrix b, NumericMatrix x0) {
  NumericMatrix x = clone(x0);
  marlib::gsstep(marlib::TRANS(), alpha, A, sigma, omega, b, x,
                 marlib::CSCMatrixT(), marlib::DenseMatrixT());
  return x;
}
