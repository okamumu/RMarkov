#ifndef MRLIB_MATRIX_H
#define MRLIB_MATRIX_H

namespace marlib {

  template <>
  struct vector_traits<double*> {
    static const double* value(const double* v) { return &v[0]; }
    static double* value(double* v) { return &v[0]; }
    static int inc(const double* v) { return 1; }
  };

  template <>
  struct vector_traits<const double*> {
    static const double* value(const double* v) { return &v[0]; }
    static int inc(const double* v) { return 1; }
  };

  template <>
  struct vector_traits<std::vector<double>> {
    static int size(const std::vector<double>& v) { return v.size(); }
    static const double* value(const std::vector<double>& v) { return v.data(); }
    static double* value(std::vector<double>& v) { return v.data(); }
    static int inc(const std::vector<double>& v) { return 1; }
  };

  class dense_matrix {
  private:
    const int m_m;
    const int m_n;
    const int m_ld;
    std::vector<double> m_value;

  public:
    dense_matrix(int m, int n)
      : m_m(m), m_n(n), m_ld(m), m_value(std::vector<double>(m*n)) {}
    dense_matrix(const dense_matrix&) = delete;
    virtual ~dense_matrix() {}
    int size() const { return m_m * m_n; }
    int nrow() const { return m_m; }
    int ncol() const { return m_n; }
    int ld() const { return m_ld; }
    const double* data() const { return m_value.data(); }
    double* data() { return m_value.data(); }
  };

  template <>
  struct vector_traits<dense_matrix> {
    static int size(const dense_matrix& v) { return v.size(); }
    static const double* value(const dense_matrix& v) { return v.data(); }
    static double* value(dense_matrix& v) { return v.data(); }
    static int inc(const dense_matrix& v) { return 1; }
  };

  template <>
  struct dense_matrix_traits<dense_matrix> {
    static int nrow(const dense_matrix& m) { return m.nrow(); }
    static int ncol(const dense_matrix& m) { return m.ncol(); }
    static const double* value(const dense_matrix& m) { return m.data(); }
    static double* value(dense_matrix& m) { return m.data(); }
    static int ld(const dense_matrix& m) { return m.ld(); }
  };

}

#endif
