#include "marlib_Rcpp.h"
using namespace Rcpp;

/**
* MatT: Type of matrix Q
* VecT1: Type of x0
* VecT2: Type of rwd
* s: the size of y = x0 * rwd
*/

// template <typename T1, typename T2, typename T3>
// void ex_ddot(const T1& x1, const T2& x2, T3& y, marlib::ArrayT, marlib::ArrayT) {
//   using traits3 = marlib::vector_traits<T3>;
//   double* valueY = traits3::value(y);
//   *valueY = marlib::ddot(x1, x2);
// }
//
// template <typename T1, typename T2, typename T3>
// void ex_ddot(const T1& x1, const T2& x2, T3& y, marlib::ArrayT, marlib::DenseMatrixT) {
//   marlib::dgemv(marlib::TRANS(), 1.0, x2, x1, 0.0, y, marlib::DenseMatrixT());
// }
//
// template <typename T1, typename T2, typename T3>
// void ex_ddot(const T1& x1, const T2& x2, T3& y, marlib::DenseMatrixT, marlib::DenseMatrixT) {
//   using traits1 = marlib::dense_matrix_traits<T1>;
//   using traits2 = marlib::dense_matrix_traits<T2>;
//   using traits3 = marlib::vector_traits<T3>;
//   const int m = traits2::ncol(x2);
//   const int k = traits2::nrow(x2);
//   const int n = traits1::ncol(x1);
//   const int ld1 = traits1::ld(x1);
//   const int ld2 = traits2::ld(x2);
//   const int ld3 = m;
//   marlib::f77dgemm('T', 'N', m, n, k, 1.0, traits2::value(x2), ld2,
//            traits1::value(x1), ld1, 0.0, traits3::value(y), ld3);
// }

template <typename TR, typename T1, typename T2, typename T3, typename T4,
          typename MatT, typename VecT1, typename VecT2>
List tran_unif_rwd(TR, const T1& Q, const T2& x0, const T2& cx0, const T3& rwd, const T4& t,
                   double ufact, double eps, int rmax, MatT, VecT1, VecT2) {
  using traits1 = marlib::vector_traits<T4>;
  const int m = traits1::size(t);
  const int s = marlib::get_size(x0, rwd, VecT1(), VecT2());
  NumericVector res_irwd(s*m);
  NumericVector res_crwd(s*m);
  T1 P = clone(Q);
  T2 x = clone(x0);
  T2 cx = clone(cx0);
  marlib::mexp_params params(ufact, eps, rmax);
  marlib::ctmc_tran_rwd(TR(), P, x, cx, rwd, t, res_irwd, res_crwd,
                        params,
                        [](marlib::mexp_params params){ stop("Time interval is too large: right = %d (rmax: %d).", params.r, params.rmax); },
                        [](marlib::mexp_params){ R_CheckUserInterrupt(); },
                        MatT(), VecT1(), VecT2());
  return List::create(
    Named("x") = x,
    Named("cx") = cx,
    Named("irwd") = res_irwd,
    Named("crwd") = res_crwd
  );
}

//////////

// [[Rcpp::export]]
List Ctran_unif_rwd_vec_vec(bool trans, S4 A, NumericVector x0, NumericVector cx0, NumericVector rwd, NumericVector t,
                        double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::DenseMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::DenseMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSRMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSRMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSCMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSCMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::COOMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::COOMatrixT(), marlib::ArrayT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Ctran_unif_rwd_vec_mat(bool trans, S4 A, NumericVector x0, NumericVector cx0, NumericMatrix rwd, NumericVector t,
                        double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::DenseMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::DenseMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSRMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSRMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSCMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSCMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::COOMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::COOMatrixT(), marlib::ArrayT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Ctran_unif_rwd_mat_mat(bool trans, S4 A, NumericMatrix x0, NumericMatrix cx0, NumericMatrix rwd, NumericVector t,
                            double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::DenseMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::DenseMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSRMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSRMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSCMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::CSCMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return tran_unif_rwd(marlib::NOTRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::COOMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return tran_unif_rwd(marlib::TRANS(), A, x0, cx0, rwd,
                           t, ufact, eps, rmax,
                           marlib::COOMatrixT(), marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}
