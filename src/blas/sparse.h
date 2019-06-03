/*
 utils for sparse matrix
*/

#pragma once

#include <cmath>

#include "traits.h"

namespace marlib {

  template <typename T1>
  int sparse_nnz(const T1& A) {
    using traits1 = dense_matrix_traits<T1>;
    const int m = traits1::nrow(A);
    const int n = traits1::ncol(A);
    const double* value = traits1::value(A);
    const int lda = traits1::ld(A);
    int nnz = 0;
    for (int j=0, cn=0; j<n; j++, cn+=lda) {
      for (int i=0; i<m; i++) {
        if (std::fpclassify(value[i+cn]) != FP_ZERO) {
          nnz++;
        }
      }
    }
    return nnz;
  }

}


