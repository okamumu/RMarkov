#pragma once

#include "traits.h"

namespace marlib {

template<typename MatT, typename T1, typename T2>
double unif(T1& Q, double ufactor, T2& tmp, MatT) {
  using traits1 = vector_traits<T2>;
  const int n = traits1::size(tmp);
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

/**
 * Compute the following:
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

/**
 * Compute the following:
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
}
