#pragma once

#include "traits.h"

namespace marlib {

inline
double& elemA(int i, int j, double* value, int ld) {
  return value[i-1+(j-1)*ld];
}

template <typename T1, typename T2>
void gth_impl(T1& A, T2& x) {
  using traits1 = vector_traits<T2>;
  using traits2 = dense_matrix_traits<T1>;
  const int n = traits1::size(x);
  double* valueX = traits1::value(x);
  const int incx = traits1::inc(x);
  // const int m = traits2::nrow(A);
  // const int n = traits2::ncol(A);
  double* valueA = traits2::value(A);
  const int lda = traits2::ld(A);

  for (int l=n; l>=2; l--) {
    // tmp = sum(A(l,1:l-1))
    double tmp = 0.0;
    for (int j=1; j<=l-1; j++) {
      tmp += elemA(l,j,valueA,lda);
    }
    for (int j=1; j<=l-1; j++) {
      for (int i=1; i<=l-1; i++) {
        if (i != j) {
          elemA(i,j,valueA,lda) += elemA(l,j,valueA,lda) * elemA(i,l,valueA,lda) / tmp;
        }
      }
    }
    for (int i=1; i<=l-1; i++) {
      elemA(i,l,valueA,lda) /= tmp;
    }
    for (int j=1; j<=l-1; j++) {
      elemA(l,j,valueA,lda) = 0.0;
    }
    elemA(l,l,valueA,lda) = -1;
  }

  double total = 0.0;
  double* tmpX = valueX;
  *tmpX = 1.0;
  total += *tmpX;
  tmpX += incx;
  for (int l=2; l<=n; l++, tmpX+=incx) {
    *tmpX = 0.0;
    double* Xi = valueX;
    for (int i=1; i<=l-1; i++, Xi+=incx) {
      *tmpX += *Xi * elemA(i,l,valueA,lda);
    }
    total += *tmpX;
  }
  dscal(1.0/total, x);
}

}
