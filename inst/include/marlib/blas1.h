#pragma once

namespace marlib {

  template<typename T1, typename T2>
  double ddot(const T1& x, const T2& y) {
    using traits1 = vector_traits<T1>;
    using traits2 = vector_traits<T2>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    const int inc2 = traits2::inc(y);
    return f77ddot(n, traits1::value(x), inc1, traits2::value(y), inc2);
  }

  template<typename T>
  double dasum(const T& x) {
    using traits1 = vector_traits<T>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    return f77dasum(n, traits1::value(x), inc1);
  }

  template<typename T>
  double dsum(const T& x) {
    using traits1 = vector_traits<T>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    return f77dsum(n, traits1::value(x), inc1);
  }

  template<typename T>
  int idamax(const T& x) {
    using traits1 = vector_traits<T>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    return f77idamax(n, traits1::value(x), inc1);
  }

  template<typename T1, typename T2>
  void dcopy(const T1& x, T2& y) {
    using traits1 = vector_traits<T1>;
    using traits2 = vector_traits<T2>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    const int inc2 = traits2::inc(y);
    return f77dcopy(n, traits1::value(x), inc1, traits2::value(y), inc2);
  }

  template<typename T>
  void dscal(double alpha, T& x) {
    using traits1 = vector_traits<T>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    f77dscal(n, alpha, traits1::value(x), inc1);
  }

  template<typename T1, typename T2>
  void daxpy(double alpha, const T1& x, T2& y) {
    using traits1 = vector_traits<T1>;
    using traits2 = vector_traits<T2>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    const int inc2 = traits2::inc(y);
    f77daxpy(n, alpha, traits1::value(x), inc1, traits2::value(y), inc2);
  }

  template<typename T>
  void dfill(T& x, double v) {
    using traits1 = vector_traits<T>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    double* p = traits1::value(x);
    for (int i=0; i<n; i++, p+=inc1) {
      *p = v;
    }
  }

  template<typename T>
  void dfill(T& x, double v, ArrayT) {
    using traits1 = vector_traits<T>;
    const int n = traits1::size(x);
    const int inc1 = traits1::inc(x);
    double* p = traits1::value(x);
    for (int i=0; i<n; i++, p+=inc1) {
      *p = v;
    }
  }


}
