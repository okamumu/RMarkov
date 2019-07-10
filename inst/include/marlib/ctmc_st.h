#ifndef MARLIB_CTMC_ST_H
#define MARLIB_CTMC_ST_H

namespace marlib {

inline
double& elemA(int i, int j, double* value, int ld) {
  return value[i-1+(j-1)*ld];
}

template <typename T1, typename T2, typename MatT>
void ctmc_st_gth(const T1& Q, T2& x, MatT) {
  const int n = nrow(Q, MatT());
  dense_matrix A(n,n);
  dcopy(Q, A, MatT(), DenseMatrixT());

  using traits1 = vector_traits<T2>;
  double* valueX = traits1::value(x);
  const int incx = traits1::inc(x);

  using traits2 = dense_matrix_traits<dense_matrix>;
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

template <typename T1, typename T2, typename Params, typename Func, typename MatT>
void ctmc_st_power(T1& P, T2& x, Params& params, const Func& callback, MatT) {
  using traits1 = vector_traits<T2>;
  int n = traits1::size(x);
  std::vector<double> tmpv(n);
  std::vector<double> prevx(n);
  unif(P, params.ufact, MatT());
  params.iter = 0;
  params.info = 1;
  while(1) {
    dcopy(x, prevx);
    for (int i=0; i<params.steps; i++) {
      dcopy(x, tmpv);
      dgemv(TRANS(), 1.0, P, tmpv, 0.0, x, MatT());
      double tmp = dasum(x);
      dscal(1.0/tmp, x);
    }
    params.rerror = drerr(prevx, x);
    params.iter += params.steps;
    if (params.rerror < params.rtol) {
      params.info = 0;
      break;
    }
    if (params.iter >= params.maxiter) {
      params.info = -1;
      break;
    }
    callback(params);
  }
}

template <typename T1, typename T2, typename Params, typename Func, typename MatT>
void ctmc_st_gs(const T1& Q, T2& x, Params& params, const Func& callback, MatT) {
  using traits1 = vector_traits<T2>;
  int n = traits1::size(x);
  std::vector<double> b(n, 0.0);
  std::vector<double> prevx(n);
  params.iter = 0;
  params.info = 1;
  while(1) {
    dcopy(x, prevx);
    for (int i=0; i<params.steps; i++) {
      gsstep(TRANS(), 1.0, Q, 0.0, 1.0, b, x, MatT(), ArrayT());
      double tmp = dasum(x);
      dscal(1.0/tmp, x);
    }
    params.rerror = drerr(prevx, x);
    params.iter += params.steps;
    if (params.rerror < params.rtol) {
      params.info = 0;
      break;
    }
    if (params.iter >= params.maxiter) {
      params.info = -1;
      break;
    }
    callback(params);
  }
}

}

#endif

