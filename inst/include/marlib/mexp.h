#ifndef MRLIB_MEXP_H
#define MRLIB_MEXP_H

namespace marlib {

template<typename MatT, typename T1, typename T2>
double unif(T1& Q, double ufactor, T2& tmp, MatT) {
  const int n = marlib::nrow(Q, MatT());
  using traits1 = vector_traits<T2>;
  double* value_tmp = traits1::value(tmp);
  const int inc = traits1::inc(tmp);
  diag_get(Q, tmp, MatT());
  double max_value = -tmp[idamax(tmp)];
  double qv = max_value * ufactor;
  dscal(1/qv, Q);
  for (int i=0; i<n; i++, value_tmp+=inc) {
    *value_tmp = 1.0 + *value_tmp/qv;
  }
  diag_set(Q, tmp, MatT());
  return qv;
}

template<typename MatT, typename T1>
double unif(T1& Q, double ufactor, MatT) {
  const int n = marlib::nrow(Q, MatT());
  std::vector<double> tmp(n);
  return unif(Q, ufactor, tmp, MatT());
}

/**
 * Compute the following: x, y are vectors
 *
 * y = exp(tr(Q) t) * x
 *
 */

template<typename MatT, typename TR, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void mexp(TR, const T1& P, const T2& poi, int right, double weight,
          const T3& x, T4& y, T5& xi, T6& tmp, MatT, ArrayT) {
  dcopy(x, xi);
  dfill(y, 0.0);
  daxpy(poi[0]/weight, xi, y);
  for (int k=1; k<=right; k++) {
    dcopy(xi, tmp);
    dgemv(TR(), 1.0, P, tmp, 0.0, xi, MatT());
    daxpy(poi[k]/weight, xi, y);
    if (fpclassify(dasum(xi)) == FP_ZERO) break;
  }
  // dscal(1.0/weight, y);
}

template<typename MatT, typename TR, typename T1, typename T2, typename T3, typename T4>
void mexp(TR, const T1& P, const T2& poi, int right, double weight,
          const T3& x, T4& y, MatT, ArrayT) {
  const int m = marlib::nrow(P, MatT());
  std::vector<double> xi(m);
  std::vector<double> tmp(m);
  mexp(TR(), P, poi, right, weight, x, y, xi, tmp, MatT(), ArrayT());
}

/**
 * Compute the following: x, y are matrices
 *
 * y = exp(tr(Q) t) * x
 *
 */

template<typename MatT, typename TR, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void mexp(TR, const T1& P, const T2& poi, int right, double weight,
          const T3& x, T4& y, T5& xi, T6& tmp, MatT, DenseMatrixT) {
  dcopy(x, xi);
  dfill(y, 0.0);
  daxpy(poi[0]/weight, xi, y);
  for (int k=1; k<=right; k++) {
    dcopy(xi, tmp);
    dgemm(TR(), NOTRANS(), 1.0, P, tmp, 0.0, xi, MatT(), DenseMatrixT());
    daxpy(poi[k]/weight, xi, y);
    if (fpclassify(dasum(xi)) == FP_ZERO) break;
  }
  // dscal(1.0/weight, y);
}

template<typename MatT, typename TR, typename T1, typename T2, typename T3, typename T4>
void mexp(TR, const T1& P, const T2& poi, int right, double weight,
          const T3& x, T4& y, MatT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T3>;
  const int m = traits1::nrow(x);
  const int n = traits1::ncol(x);
  dense_matrix xi(m,n);
  dense_matrix tmp(m,n);
  mexp(TR(), P, poi, right, weight, x, y, xi, tmp, MatT(), DenseMatrixT());
}

/**
 * Compute the following: x, cy, y are vector
 *
 * cy = cy + int_0^t exp(tr(Q) s) * x ds
 * y = y + exp(tr(Q) s) * x ds
 *
 */

template<typename MatT, typename TR,
         typename T1, typename T2, typename T3, typename T4,
         typename T5, typename T6, typename T7, typename T8>
void mexpint(TR, const T1& P, const T2& poi, const T3& cpoi, int right, double weight, double qv_weight,
             const T4& x, T5& y, T6& cy, T7& xi, T8& tmp, MatT, ArrayT) {
  dcopy(x, xi);
  dfill(y, 0.0);
  daxpy(poi[0]/weight, xi, y);
  daxpy(cpoi[0]/qv_weight, xi, cy);
  for (int k=1; k<=right; k++) {
    dcopy(xi, tmp);
    dgemv(TR(), 1.0, P, tmp, 0.0, xi, MatT());
    daxpy(poi[k]/weight, xi, y);
    daxpy(cpoi[k]/qv_weight, xi, cy);
    if (fpclassify(dasum(xi)) == FP_ZERO) break;
  }
}

template<typename MatT, typename TR,
         typename T1, typename T2, typename T3, typename T4,
         typename T5, typename T6>
void mexpint(TR, const T1& P, const T2& poi, const T3& cpoi, int right, double weight, double qv_weight,
             const T4& x, T5& y, T6& cy, MatT, ArrayT) {
  const int m = nrow(P, MatT());
  std::vector<double> xi(m);
  std::vector<double> tmp(m);
  mexpint(TR(), P, poi, cpoi, right, weight, qv_weight,
          x, y, cy, xi, tmp, MatT(), ArrayT());
}

/**
 * Compute the following: x, cy, y are matrix
 *
 * cy = cy + int_0^t exp(tr(Q) s) * x ds
 * y = y + exp(tr(Q) s) * x ds
 *
 */

template<typename MatT, typename TR,
         typename T1, typename T2, typename T3, typename T4,
         typename T5, typename T6, typename T7, typename T8>
void mexpint(TR, const T1& P, const T2& poi, const T3& cpoi, int right, double weight, double qv_weight,
             const T4& x, T5& y, T6& cy, T7& xi, T8& tmp, MatT, DenseMatrixT) {
  dcopy(x, xi);
  dfill(y, 0.0);
  daxpy(poi[0]/weight, xi, y);
  daxpy(cpoi[0]/qv_weight, xi, cy);
  for (int k=1; k<=right; k++) {
    dcopy(xi, tmp);
    dgemm(TR(), NOTRANS(), 1.0, P, tmp, 0.0, xi, MatT(), DenseMatrixT());
    daxpy(poi[k]/weight, xi, y);
    daxpy(cpoi[k]/qv_weight, xi, cy);
    if (fpclassify(dasum(xi)) == FP_ZERO) break;
  }
}

template<typename MatT, typename TR,
         typename T1, typename T2, typename T3, typename T4,
         typename T5, typename T6>
void mexpint(TR, const T1& P, const T2& poi, const T3& cpoi, int right, double weight, double qv_weight,
             const T4& x, T5& y, T6& cy, MatT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T4>;
  const int m = traits1::nrow(x);
  const int n = traits1::ncol(x);
  dense_matrix xi(m,n);
  dense_matrix tmp(m,n);
  mexpint(TR(), P, poi, cpoi, right, weight, qv_weight,
          x, y, cy, xi, tmp, MatT(), DenseMatrixT());
}

// template<typename TR>
// struct not_{};
//
// template<>
// struct not_<trans>{
//   using Type = notrans;
// };
//
// template<>
// struct not_<notrans>{
//   using Type = trans;
// };
//
// //
// template<typename MatT, typename TR,
//          typename T1, typename T2, typename T3, typename T4,
//          typename T5, typename T6, typename T7, typename T8>
// void mexp_conv(MatT, TR, const T1& P, double qv,
//                const T2& poi, int right, double weight,
//                const T3& x, const T4& y, T5& z, T6& H, T7& xi, T8& vc) {
//   using notTR = not_<TR>;
//   for (int i=0; i<=right+1; i++) {
//     dfill(vc[i], 0.0);
//   }
//
//   // dfill(n, vc[right+1], 0.0);
//   daxpy(poi[right+1], y, vc[right+1]);
//   for (int l=right; l>=1; l--) {
//     dgemv(MatT(), typename notTR::Type(), 1.0, P, vc[l+1], 0.0, vc[l]);
//     daxpy(poi[l], y, vc[l]);
//   }
//
//   dcopy(x, xi);
//   dfill(z, 0.0);
//   daxpy(poi[0], xi, z);
//   // dfill(H, 0.0);
//   dger(MatT(), 1.0/qv/weight, xi, vc[1], H);
//   for (int l=1; l<=right; l++) {
//     dgemv(MatT(), TR(), 1.0, P, xi, 0.0, xi);
//     daxpy(poi[l], xi, z);
//     dger(MatT(), 1.0/qv/weight, xi, vc[l+1], H);
//   }
//   dscal(1.0/weight, z);
//   // dscal(1.0/qv/weight, H);
// }
//
// }

//// high level

template <typename TR, typename T1, typename T2, typename T3,
          typename Params, typename Func, typename MatT, typename VecT>
void mexp_func(TR, T1& P, const T2& x, T3& y,
               double t, Params& params, const Func& callback, MatT, VecT) {
  double qv = unif(P, params.ufact, MatT());
  params.r = poi::rightbound(qv*t, params.eps);
  if (params.r > params.rmax) {
    callback(params); //stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  std::vector<double> prob(params.r+1);
  double weight = poi::pmf(qv*t, 0, params.r, prob);
  mexp(TR(), P, prob, params.r, weight, x, y, MatT(), VecT());
}

template <typename TR, typename T1, typename T2, typename T3, typename T4,
          typename Params, typename Func, typename MatT, typename VecT>
void mexpint_func(TR, T1& P, const T2& x, T3& y, T4& cy,
                  double t, Params& params, const Func& callback, MatT, VecT) {
  double qv = unif(P, params.ufact, MatT());
  params.r = poi::rightbound(qv*t, params.eps) + 1;
  if (params.r > params.rmax) {
    callback(params); //stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  std::vector<double> prob(params.r+1);
  std::vector<double> cprob(params.r+1);
  double weight = poi::cpmf(qv*t, 0, params.r, prob, cprob);
  mexpint(TR(), P, prob, cprob, params.r, weight, qv*weight, x, y, cy, MatT(), VecT());
}

}

#endif
