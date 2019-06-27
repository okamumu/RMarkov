#include "marlib_Rcpp.h"
using namespace Rcpp;

// general

template <typename T1, typename MatT>
List C_unif(int n, T1 Q, double ufactor) {
  NumericVector tmp(n);
  T1 P = clone(Q);
  double qv = marlib::unif(P, ufactor, tmp, MatT());
  return List::create(Named("P") = P, Named("qv")=qv);
}

template <typename TR, typename T1, typename MatT>
NumericVector C_mexp_vec(T1 Q, double t, NumericVector x, double eps, double ufactor) {
  const int n = x.size();
  NumericVector xi(n);
  NumericVector tmp(n);
  NumericVector y(n);
  T1 P = clone(Q);
  double qv = marlib::unif(P, ufactor, tmp, MatT());
  int r = marlib::poi::rightbound(qv*t, eps);
  NumericVector prob(r+1);
  double weight = marlib::poi::pmf(qv*t, 0, r, prob);
  marlib::mexp(TR(), P, prob, r, weight, x, y, xi, tmp, MatT(), marlib::ArrayT());
  return y;
}

template <typename TR, typename T1, typename MatT>
NumericMatrix C_mexp_mat(T1 Q, double t, NumericMatrix x, double eps, double ufactor) {
  const int m = x.nrow();
  const int n = x.ncol();
  NumericMatrix xi(m,n);
  NumericMatrix tmp(m,n);
  NumericMatrix y(m,n);
  T1 P = clone(Q);
  double qv = marlib::unif(P, ufactor, tmp, MatT());
  int r = marlib::poi::rightbound(qv*t, eps);
  NumericVector prob(r+1);
  double weight = marlib::poi::pmf(qv*t, 0, r, prob);
  marlib::mexp(TR(), P, prob, r, weight, x, y, xi, tmp, MatT(), marlib::DenseMatrixT());
  return y;
}

template <typename TR, typename T1, typename MatT>
List C_mexpint_vec(T1 Q, double t, NumericVector x, NumericVector cy0, double eps, double ufactor) {
  const int n = x.size();
  NumericVector xi(n);
  NumericVector tmp(n);
  NumericVector y(n);
  T1 P = clone(Q);
  NumericVector cy = clone(cy0);
  double qv = marlib::unif(P, ufactor, tmp, MatT());
  int r = marlib::poi::rightbound(qv*t, eps);
  NumericVector prob(r+1);
  NumericVector cprob(r+1);
  double weight = marlib::poi::cpmf(qv*t, 0, r, prob, cprob);
  marlib::mexpint(TR(), P, prob, cprob, r, weight, qv*weight, x, y, cy, xi, tmp, MatT(), marlib::ArrayT());
  return List::create(Named("y")=y, Named("cy")=cy);
}

template <typename TR, typename T1, typename MatT>
List C_mexpint_mat(T1 Q, double t, NumericMatrix x, NumericMatrix cy0, double eps, double ufactor) {
  const int m = x.nrow();
  const int n = x.ncol();
  NumericMatrix xi(m,n);
  NumericMatrix tmp(m,n);
  NumericMatrix y(m,n);
  T1 P = clone(Q);
  NumericMatrix cy = clone(cy0);
  double qv = marlib::unif(P, ufactor, tmp, MatT());
  int r = marlib::poi::rightbound(qv*t, eps);
  NumericVector prob(r+1);
  NumericVector cprob(r+1);
  double weight = marlib::poi::cpmf(qv*t, 0, r, prob, cprob);
  marlib::mexpint(TR(), P, prob, cprob, r, weight, qv*weight, x, y, cy, xi, tmp, MatT(), marlib::DenseMatrixT());
  return List::create(Named("y")=y, Named("cy")=cy);
}

//' @export
// [[Rcpp::export]]
List C_unif_dense(NumericMatrix Q, double ufactor) {
  const int n = Q.nrow();
  return C_unif<NumericMatrix, marlib::DenseMatrixT>(n, Q, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_unif_csr(S4 Q, double ufactor) {
  const int n = Rcpp::as<Rcpp::IntegerVector>(Q.slot("Dim"))[1];
  return C_unif<S4, marlib::CSRMatrixT>(n, Q, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_unif_csc(S4 Q, double ufactor) {
  const int n = Rcpp::as<Rcpp::IntegerVector>(Q.slot("Dim"))[1];
  return C_unif<S4, marlib::CSCMatrixT>(n, Q, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_unif_coo(S4 Q, double ufactor) {
  const int n = Rcpp::as<Rcpp::IntegerVector>(Q.slot("Dim"))[1];
  return C_unif<S4, marlib::COOMatrixT>(n, Q, ufactor);
}

///// mexp

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecT_dense(NumericMatrix Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::TRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecN_dense(NumericMatrix Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::NOTRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecT_csr(S4 Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::TRANS, S4, marlib::CSRMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecN_csr(S4 Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::NOTRANS, S4, marlib::CSRMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecT_csc(S4 Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::TRANS, S4, marlib::CSCMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecN_csc(S4 Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::NOTRANS, S4, marlib::CSCMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecT_coo(S4 Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::TRANS, S4, marlib::COOMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericVector C_mexp_vecN_coo(S4 Q, double t, NumericVector x, double eps, double ufactor) {
  return C_mexp_vec<marlib::NOTRANS, S4, marlib::COOMatrixT>(Q, t, x, eps, ufactor);
}

///// mexp mat

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matN_dense(NumericMatrix Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::NOTRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matT_dense(NumericMatrix Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::TRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matN_csr(S4 Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::NOTRANS, S4, marlib::CSRMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matT_csr(S4 Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::TRANS, S4, marlib::CSRMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matN_csc(S4 Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::NOTRANS, S4, marlib::CSCMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matT_csc(S4 Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::TRANS, S4, marlib::CSCMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matN_coo(S4 Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::NOTRANS, S4, marlib::COOMatrixT>(Q, t, x, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_mexp_matT_coo(S4 Q, double t, NumericMatrix x, double eps, double ufactor) {
  return C_mexp_mat<marlib::TRANS, S4, marlib::COOMatrixT>(Q, t, x, eps, ufactor);
}

/// mexpint

//' @export
// [[Rcpp::export]]
List C_mexpint_vecT_dense(NumericMatrix Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::TRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecN_dense(NumericMatrix Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::NOTRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecT_csr(S4 Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::TRANS, S4, marlib::CSRMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecN_csr(S4 Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::NOTRANS, S4, marlib::CSRMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecT_csc(S4 Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::TRANS, S4, marlib::CSCMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecN_csc(S4 Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::NOTRANS, S4, marlib::CSCMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecT_coo(S4 Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::TRANS, S4, marlib::COOMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_vecN_coo(S4 Q, double t, NumericVector x, NumericVector cy, double eps, double ufactor) {
  return C_mexpint_vec<marlib::NOTRANS, S4, marlib::COOMatrixT>(Q, t, x, cy, eps, ufactor);
}

/// mexpint mat

//' @export
// [[Rcpp::export]]
List C_mexpint_matT_dense(NumericMatrix Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::TRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matN_dense(NumericMatrix Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::NOTRANS, NumericMatrix, marlib::DenseMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matT_csr(S4 Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::TRANS, S4, marlib::CSRMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matN_csr(S4 Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::NOTRANS, S4, marlib::CSRMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matT_csc(S4 Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::TRANS, S4, marlib::CSCMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matN_csc(S4 Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::NOTRANS, S4, marlib::CSCMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matT_coo(S4 Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::TRANS, S4, marlib::COOMatrixT>(Q, t, x, cy, eps, ufactor);
}

//' @export
// [[Rcpp::export]]
List C_mexpint_matN_coo(S4 Q, double t, NumericMatrix x, NumericMatrix cy, double eps, double ufactor) {
  return C_mexpint_mat<marlib::NOTRANS, S4, marlib::COOMatrixT>(Q, t, x, cy, eps, ufactor);
}
