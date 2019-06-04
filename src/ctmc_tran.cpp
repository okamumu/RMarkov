#include <Rcpp.h>
#include <marlib.h>

using namespace Rcpp;

template <typename TR, typename T1, typename T2, typename T3, typename MatT, typename VecT>
List tran_unif(TR, const T1& Q, const T2& x0, const T3& t,
               double ufact, double eps, int rmax, MatT, VecT) {
  using traits1 = marlib::vector_traits<T3>;
  using traits2 = marlib::vector_traits<T2>;
  const int m = traits1::size(t);
  const double* valueT = traits1::value(t);
  const int inct = traits1::inc(t);
  const int s = traits2::size(x0);
  NumericVector res_x(s*m);
  NumericVector res_cx(s*m);

  T1 P = clone(Q);
  T2 xi = clone(x0);
  T2 tmp = clone(x0);
  T2 x = clone(x0);
  T2 cx = clone(x0);
  double qv = marlib::unif(P, ufact, tmp, MatT());
  double maxt = t[marlib::idamax(t)];
  int r = marlib::poi::rightbound(qv*maxt, eps) + 1;
  if (r > rmax) {
    stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  NumericVector prob(r+1);
  NumericVector cprob(r+1);
  marlib::dfill(cx, 0.0);
  double* ptr_x = &res_x[0];
  double* ptr_cx = &res_cx[0];
  for (int i=0; i<m; i++, valueT+=inct, ptr_x+=s, ptr_cx+=s) {
    double t = *valueT;
    int r = marlib::poi::rightbound(qv*t, eps) + 1;
    double weight = marlib::poi::cpmf(qv*t, 0, r, prob, cprob);
    marlib::mexpint(TR(), P, prob, cprob, r, weight, qv*weight, x, x, cx, xi, tmp, MatT(), VecT());
    marlib::dcopy(x, ptr_x);
    marlib::dcopy(cx, ptr_cx);
    R_CheckUserInterrupt();
  }
  return List::create(
    Named("x") = res_x,
    Named("cx") = res_cx
  );
}

//////////

// [[Rcpp::export]]
List Ctran_unif_vec(bool trans, S4 A, NumericVector x0, NumericVector t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Ctran_unif_mat(bool trans, S4 A, NumericMatrix x0, NumericVector t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}
