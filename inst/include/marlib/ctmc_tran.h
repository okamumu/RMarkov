#ifndef MRLIB_CTMC_TRAN_H
#define MRLIB_CTMC_TRAN_H

namespace marlib {

template <typename TR, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
          typename Params, typename Func1, typename Func2,
          typename MatT, typename VecT>
void ctmc_tran(TR, T1& P, T2& x, T3& cx, const T4& t, T5& res_x, T6& res_cx,
               Params& params, const Func1& callback1, const Func2& callback2, MatT, VecT) {
  using traits1 = vector_traits<T4>;
  using traits2 = vector_traits<T2>;
  using traits3 = vector_traits<T5>;
  using traits4 = vector_traits<T6>;
  const int m = traits1::size(t);
  const double* valueT = traits1::value(t);
  const int inct = traits1::inc(t);
  const int s = traits2::size(x);
  double* ptr_x = traits3::value(res_x);
  double* ptr_cx = traits4::value(res_cx);

  // dcopy(x0, x);
  // dfill(cx, 0.0);
  double qv = unif(P, params.ufact, MatT());
  double maxt = dmax(t);
  params.r = poi::rightbound(qv*maxt, params.eps) + 1;
  if (params.r > params.rmax) {
    callback1(params); //stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  std::vector<double> prob(params.r+1);
  std::vector<double> cprob(params.r+1);
  for (int i=0; i<m; i++, valueT+=inct, ptr_x+=s, ptr_cx+=s) {
    double t = *valueT;
    int r = poi::rightbound(qv*t, params.eps) + 1;
    double weight = poi::cpmf(qv*t, 0, r, prob, cprob);
    mexpint(TR(), P, prob, cprob, r, weight, qv*weight, x, x, cx, MatT(), VecT());
    dcopy(x, ptr_x);
    dcopy(cx, ptr_cx);
    callback2(params); // R_CheckUserInterrupt();
  }
}

/**
 * MatT: Type of matrix Q
 * VecT1: Type of x0
 * VecT2: Type of rwd
 * s: the size of y = x0 * rwd
 */

template <typename T1, typename T2>
int get_size(const T1& x1, const T2& x2, ArrayT, ArrayT) {
  return 1;
}

template <typename T1, typename T2>
int get_size(const T1& x1, const T2& x2, ArrayT, DenseMatrixT) {
  using traits2 = dense_matrix_traits<T2>;
  return traits2::ncol(x2);
}

template <typename T1, typename T2>
int get_size(const T1& x1, const T2& x2, DenseMatrixT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  return traits1::ncol(x1) * traits2::ncol(x2);
}

template <typename T1, typename T2, typename T3>
void ex_ddot(const T1& x1, const T2& x2, T3& y, ArrayT, ArrayT) {
  using traits3 = vector_traits<T3>;
  double* valueY = traits3::value(y);
  *valueY = ddot(x1, x2);
}

template <typename T1, typename T2, typename T3>
void ex_ddot(const T1& x1, const T2& x2, T3& y, ArrayT, DenseMatrixT) {
  dgemv(TRANS(), 1.0, x2, x1, 0.0, y, DenseMatrixT());
}

template <typename T1, typename T2, typename T3>
void ex_ddot(const T1& x1, const T2& x2, T3& y, DenseMatrixT, DenseMatrixT) {
  using traits1 = dense_matrix_traits<T1>;
  using traits2 = dense_matrix_traits<T2>;
  using traits3 = vector_traits<T3>;
  const int m = traits2::ncol(x2);
  const int k = traits2::nrow(x2);
  const int n = traits1::ncol(x1);
  const int ld1 = traits1::ld(x1);
  const int ld2 = traits2::ld(x2);
  const int ld3 = m;
  f77dgemm('T', 'N', m, n, k, 1.0, traits2::value(x2), ld2,
                   traits1::value(x1), ld1, 0.0, traits3::value(y), ld3);
}

template <typename TR, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6, typename T7,
          typename Params, typename Func1, typename Func2,
          typename MatT, typename VecT1, typename VecT2>
void ctmc_tran_rwd(TR, T1& P, T2& x, T3& cx, const T4& rwd, const T5& t,
                   T6& res_irwd, T7& res_crwd,
                   Params& params, const Func1& callback1, const Func2& callback2,
                   MatT, VecT1, VecT2) {
  using traits1 = vector_traits<T5>;
  using traits3 = vector_traits<T6>;
  using traits4 = vector_traits<T7>;
  const int m = traits1::size(t);
  const double* valueT = traits1::value(t);
  const int inct = traits1::inc(t);
  const int s = get_size(x, rwd, VecT1(), VecT2());
  double* ptr_x = traits3::value(res_irwd);
  double* ptr_c = traits4::value(res_crwd);

  double qv = unif(P, params.ufact, MatT());
  double maxt = dmax(t);
  params.r = poi::rightbound(qv*maxt, params.eps) + 1;
  if (params.r > params.rmax) {
    callback1(params); //stop("Time interval is too large: right = %d (rmax: %d).", r, rmax);
  }
  std::vector<double> prob(params.r+1);
  std::vector<double> cprob(params.r+1);
  for (int i=0; i<m; i++, valueT+=inct, ptr_x+=s, ptr_c+=s) {
    double t = *valueT;
    int r = poi::rightbound(qv*t, params.eps) + 1;
    double weight = poi::cpmf(qv*t, 0, r, prob, cprob);
    mexpint(TR(), P, prob, cprob, r, weight, qv*weight, x, x, cx, MatT(), VecT1());
    ex_ddot(x, rwd, ptr_x, VecT1(), VecT2());
    ex_ddot(cx, rwd, ptr_c, VecT1(), VecT2());
    callback2(params); // R_CheckUserInterrupt();
  }
}

}


// template <typename TR, typename T1, typename T2, typename T3, typename T4,
//           typename MatT, typename VecT1, typename VecT2>
// List tran_unif_rwd(TR, const T1& Q, const T2& x0, const T3& rwd, int s, const T4& t,
//                    double ufact, double eps, int rmax, MatT, VecT1, VecT2) {
//   using traits1 = vector_traits<T4>;
//   const int m = traits1::size(t);
//   const double* valueT = traits1::value(t);
//   const int inct = traits1::inc(t);
//   NumericVector res_irwd(s*m);
//   NumericVector res_crwd(s*m);
//
//   T1 P = clone(Q);
//   T2 xi = clone(x0);
//   T2 tmp = clone(x0);
//   T2 x = clone(x0);
//   T2 cx = clone(x0);
//
#endif

