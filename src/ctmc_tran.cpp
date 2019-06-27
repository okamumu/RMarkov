#include "marlib_Rcpp.h"
using namespace Rcpp;

template <typename TR, typename T1, typename T2, typename T3, typename MatT, typename VecT>
List tran_unif(TR, const T1& Q, const T2& x0, const T2& cx0, const T3& t,
               double ufact, double eps, int rmax, MatT, VecT) {
  using traits1 = marlib::vector_traits<T3>;
  using traits2 = marlib::vector_traits<T2>;
  const int m = traits1::size(t);
  const int s = traits2::size(x0);
  NumericVector res_x(s*m);
  NumericVector res_cx(s*m);
  T1 P = clone(Q);
  T2 x = clone(x0);
  T2 cx = clone(cx0);
  marlib::mexp_params params(ufact, eps, rmax);
  marlib::ctmc_tran(TR(), P, x, cx, t, res_x, res_cx,
                    params,
                    [](marlib::mexp_params params){ stop("Time interval is too large: right = %d (rmax: %d).", params.r, params.rmax); },
                    [](marlib::mexp_params){ R_CheckUserInterrupt(); },
                    MatT(), VecT());
  return List::create(
    Named("x") = res_x,
    Named("cx") = res_cx
  );
}

//////////

// [[Rcpp::export]]
List Ctran_unif_vec(bool trans, S4 A, NumericVector x0, NumericVector cx0,
                    NumericVector t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Ctran_unif_mat(bool trans, S4 A, NumericMatrix x0, NumericMatrix cx0, NumericVector t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif(marlib::NOTRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif(marlib::TRANS(), A, x0, cx0, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}
