#ifndef MARLIB_CTMC_STSEN_H
#define MARLIB_CTMC_STSEN_H

namespace marlib {

template <typename T1, typename T2, typename T3, typename T4,
          typename Params, typename Func, typename MatT>
void ctmc_stsen_gs(const T1& Q, T2& x, const T3& b, const T4& pis,
                   Params& params, const Func& callback, MatT) {
  using traits1 = vector_traits<T2>;
  int n = traits1::size(x);
  std::vector<double> prevx(n);
  params.iter = 0;
  params.info = 1;
  while(1) {
    dcopy(x, prevx);
    for (int i=0; i<params.steps; i++) {
      gsstep(TRANS(), -1.0, Q, 0.0, 1.0, b, x, MatT(), ArrayT());
      double tmp = dsum(x);
      daxpy(-tmp, pis, x);
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

