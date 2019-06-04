#include <Rcpp.h>
#include <marlib.h>

using namespace Rcpp;

template <typename TR, typename T1, typename T2, typename MatT, typename VecT>
T2 mexp_unif(TR, const T1& A, T2 x, double t, double ufact, double eps, int rmax, MatT, VecT) {
  T1 P = clone(A);
  T2 xi = clone(x);
  T2 tmp = clone(x);
  T2 y = clone(x);
  double qv = marlib::unif(P, ufact, tmp, MatT());
  int r = marlib::poi::rightbound(qv*t, eps);
  if (r > rmax) {
    stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  NumericVector prob(r+1);
  double weight = marlib::poi::pmf(qv*t, 0, r, prob);
  marlib::mexp(TR(), P, prob, r, weight, x, y, xi, tmp, MatT(), VecT());
  return y;
}

template <typename TR, typename T1, typename T2, typename MatT, typename VecT>
List mexpint_unif(TR, const T1& A, T2 x, double t, double ufact, double eps, int rmax, MatT, VecT) {
  T1 P = clone(A);
  T2 xi = clone(x);
  T2 tmp = clone(x);
  T2 y = clone(x);
  T2 cy = clone(x);
  double qv = marlib::unif(P, ufact, tmp, MatT());
  int r = marlib::poi::rightbound(qv*t, eps) + 1;
  if (r > rmax) {
    stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  NumericVector prob(r+1);
  NumericVector cprob(r+1);
  double weight = marlib::poi::cpmf(qv*t, 0, r, prob, cprob);
  marlib::mexpint(TR(), P, prob, cprob, r, weight, qv*weight, x, y, cy, xi, tmp, MatT(), VecT());
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
List Cmexpint_unif_vec(bool trans, S4 A, NumericVector x, double t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::ArrayT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::ArrayT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}

// [[Rcpp::export]]
List Cmexpint_unif_mat(bool trans, S4 A, NumericMatrix x, double t, double ufact, double eps, int rmax) {
  std::string classname = as<std::string>(A.attr("class"));
  if (classname == "dgeMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::DenseMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgRMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSRMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgCMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::CSCMatrixT(), marlib::DenseMatrixT());
    }
  } else if (classname == "dgTMatrix") {
    if (trans == false) {
      return mexpint_unif(marlib::NOTRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    } else {
      return mexpint_unif(marlib::TRANS(), A, x, t, ufact, eps, rmax, marlib::COOMatrixT(), marlib::DenseMatrixT());
    }
  } else {
    stop("Invalid class: %s", classname);
  }
}
