#include "marlib_Rcpp.h"
using namespace Rcpp;

template <typename TR, typename T1, typename T2, typename T3, typename T4,
          typename MatT, typename VecT>
T2 mexp_mix_unif(TR, const T1& A, const T2& x, const T3& w, const T4& t,
                   double ufact, double eps, int rmax, MatT, VecT) {
  T1 P = clone(A);
  T2 v = clone(x);
  T2 y = clone(x);
  marlib::mexp_params params(ufact, eps, rmax);
  marlib::mexp_mix(TR(), P, v, y, w, t, params,
                   [](marlib::mexp_params params){ stop("Time interval is too large: right = %d (rmax: %d).", params.r, params.rmax); },
                   [](marlib::mexp_params){ R_CheckUserInterrupt(); },
                   MatT(), VecT());
  return y;
}

template <typename TR, typename T1, typename T2, typename T3, typename T4,
          typename MatT, typename VecT>
List mexpint_mix_unif(TR, const T1& A, const T2& x, const T3& w, const T4& t,
                 double ufact, double eps, int rmax, MatT, VecT) {
  T1 P = clone(A);
  T2 v = clone(x);
  T2 y = clone(x);
  T2 cv = clone(x);
  T2 cy = clone(x);
  marlib::mexp_params params(ufact, eps, rmax);
  marlib::mexpint_mix(TR(), P, v, cv, y, cy, w, t, params,
                      [](marlib::mexp_params params){ stop("Time interval is too large: right = %d (rmax: %d).", params.r, params.rmax); },
                      [](marlib::mexp_params){ R_CheckUserInterrupt(); },
                      MatT(), VecT());
  return List::create(Named("y")=y, Named("cy")=cy);
}

///////////////////////////

// [[Rcpp::export]]
NumericVector Cmexp_mix_unif_vec(bool trans, S4 A, NumericVector x,
                                 NumericVector w, NumericVector t,
                                 double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
NumericMatrix Cmexp_mix_unif_mat(bool trans, S4 A, NumericMatrix x,
                                 NumericVector w, NumericVector t,
                                 double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexp_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexp_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

///////////////////////////

// [[Rcpp::export]]
List Cmexpint_mix_unif_vec(bool trans, S4 A, NumericVector x,
                           NumericVector w, NumericVector t,
                           double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Cmexpint_mix_unif_mat(bool trans, S4 A, NumericMatrix x,
                           NumericVector w, NumericVector t,
                           double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexpint_mix_unif(marlib::NOTRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_mix_unif(marlib::TRANS(), A, x, w, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}
