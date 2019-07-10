#include "marlib_Rcpp.h"
using namespace Rcpp;

template <typename TR, typename T1, typename T2, typename MatT, typename VecT>
T2 mexp_unif(TR, const T1& A, T2 x, double t, double ufact, double eps, int rmax, MatT, VecT) {
  T1 P = clone(A);
  T2 y = clone(x);
  marlib::marlib_params params;
  params.ufact = ufact;
  params.eps = eps;
  params.rmax = rmax;
  marlib::mexp_func(TR(), P, x, y, t, params,
            [](const marlib::marlib_params& params){
              stop("Time interval is too large: right = %d (rmax: %d).", params.r, params.rmax);
              },
            MatT(), VecT());
  return y;
}

template <typename TR, typename T1, typename T2, typename T3, typename MatT, typename VecT>
List mexpint_unif(TR, const T1& A, const T2& x, const T3& cx,
                  double t, double ufact, double eps, int rmax, MatT, VecT) {
  T1 P = clone(A);
  T2 y = clone(x);
  T2 cy = clone(cx);
  marlib::marlib_params params;
  params.ufact = ufact;
  params.eps = eps;
  params.rmax = rmax;
  marlib::mexpint_func(TR(), P, x, y, cy, t, params,
                      [](const marlib::marlib_params& params){
                        stop("Time interval is too large: right = %d (rmax: %d).", params.r, params.rmax);
                      },
                      MatT(), VecT());
  return List::create(Named("y")=y, Named("cy")=cy);
}

///////////////////////////

// [[Rcpp::export]]
NumericVector Cmexp_unif_vec(bool trans, S4 A, NumericVector x, double t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
NumericMatrix Cmexp_unif_mat(bool trans, S4 A, NumericMatrix x, double t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexp_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

////////////

// [[Rcpp::export]]
List Cmexpint_unif_vec(bool trans, S4 A, NumericVector x, NumericVector cx, double t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Cmexpint_unif_mat(bool trans, S4 A, NumericMatrix x, NumericMatrix cx, double t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, cx, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, cx, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}
