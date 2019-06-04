#include <Rcpp.h>
#include <marlib.h>

using namespace Rcpp;

template <typename TR, typename T1, typename T2, typename T3, typename T4,
          typename MatT, typename VecT>
T2 mexp_mix_unif(TR, const T1& A, const T2& x, const T3& w, const T4& t,
                   double ufact, double eps, int rmax, MatT, VecT) {
  using traits1 = marlib::vector_traits<T3>;
  using traits2 = marlib::vector_traits<T4>;
  const int m = traits1::size(w);
  const double* valueW = traits1::value(w);
  const double* valueT = traits2::value(t);
  const int incw = traits1::inc(w);
  const int inct = traits2::inc(t);

  T1 P = clone(A);
  T2 xi = clone(x);
  T2 tmp = clone(x);
  T2 v = clone(x);
  T2 y = clone(x);

  double qv = marlib::unif(P, ufact, tmp, MatT());
  double maxt = t[marlib::idamax(t)];
  int r = marlib::poi::rightbound(qv*maxt, eps);
  if (r > rmax) {
    stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  NumericVector prob(r+1);
  marlib::dfill(y, 0.0);
  for (int i=0; i<m; i++, valueW+=incw, valueT+=inct) {
    double w = *valueW;
    double t = *valueT;
    int r = marlib::poi::rightbound(qv*t, eps);
    double weight = marlib::poi::pmf(qv*t, 0, r, prob);
    marlib::mexp(TR(), P, prob, r, weight, v, v, xi, tmp, MatT(), VecT());
    marlib::daxpy(w, v, y);
    R_CheckUserInterrupt();
  }
  return y;
}

template <typename TR, typename T1, typename T2, typename T3, typename T4,
          typename MatT, typename VecT>
List mexpint_mix_unif(TR, const T1& A, const T2& x, const T3& w, const T4& t,
                 double ufact, double eps, int rmax, MatT, VecT) {
  using traits1 = marlib::vector_traits<T3>;
  using traits2 = marlib::vector_traits<T4>;
  const int m = traits1::size(w);
  const double* valueW = traits1::value(w);
  const double* valueT = traits2::value(t);
  const int incw = traits1::inc(w);
  const int inct = traits2::inc(t);

  T1 P = clone(A);
  T2 xi = clone(x);
  T2 tmp = clone(x);
  T2 v = clone(x);
  T2 y = clone(x);
  T2 cv = clone(x);
  T2 cy = clone(x);

  double qv = marlib::unif(P, ufact, tmp, MatT());
  double maxt = t[marlib::idamax(t)];
  int r = marlib::poi::rightbound(qv*maxt, eps) + 1;
  if (r > rmax) {
    stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  NumericVector prob(r+1);
  NumericVector cprob(r+1);
  marlib::dfill(y, 0.0);
  marlib::dfill(cv, 0.0);
  marlib::dfill(cy, 0.0);
  for (int i=0; i<m; i++, valueW+=incw, valueT+=inct) {
    double w = *valueW;
    double t = *valueT;
    int r = marlib::poi::rightbound(qv*t, eps) + 1;
    double weight = marlib::poi::cpmf(qv*t, 0, r, prob, cprob);
    marlib::mexpint(TR(), P, prob, cprob, r, weight, qv*weight, v, v, cv, xi, tmp, MatT(), VecT());
    marlib::daxpy(w, v, y);
    marlib::daxpy(w, cv, cy);
    R_CheckUserInterrupt();
  }
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
