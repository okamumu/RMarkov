#pragma once

namespace marlib {

  class dense_matrix {
  private:
    const int m_m;
    const int m_n;
    const int m_ld;

    double* m_value;

  public:
    dense_matrix(int m, int n, int ld, double* value)
      : m_m(m), m_n(n), m_ld(ld), m_value(value) {}
    dense_matrix(const dense_matrix& x)
      : m_m(x.m_m), m_n(x.m_n), m_ld(x.m_ld), m_value(x.m_value) {}
    virtual ~dense_matrix() {}
    int size() const { return m_m * m_n; }
    int nrow() const { return m_m; }
    int ncol() const { return m_n; }
    int ld() const { return m_ld; }

    double& operator[](int i) { return m_value[i]; }
    const double& operator[](int i) const { return m_value[i]; }
  };

  class const_dense_matrix {
  private:
    const int m_m;
    const int m_n;
    const int m_ld;

    const double* m_value;

  public:
    const_dense_matrix(int m, int n, int ld, const double* value)
      : m_m(m), m_n(n), m_ld(ld), m_value(value) {}
    const_dense_matrix(const const_dense_matrix& x)
      : m_m(x.m_m), m_n(x.m_n), m_ld(x.m_ld), m_value(x.m_value) {}
    virtual ~const_dense_matrix() {}
    int size() const { return m_m * m_n; }
    int nrow() const { return m_m; }
    int ncol() const { return m_n; }
    int ld() const { return m_ld; }

    const double& operator[](int i) const { return m_value[i]; }
  };

  class csr_matrix {
  private:
    const int m_m;
    const int m_n;
    const int m_nnz;
    const int m_base;

    double* m_value;
    const int* m_rowptr;
    const int* m_colind;

  public:
    csr_matrix(int m, int n, int nnz, int base, double* value, const int* rowptr, const int* colind)
    : m_m(m), m_n(n), m_nnz(nnz), m_base(base),
      m_value(value), m_rowptr(rowptr), m_colind(colind) {}
    csr_matrix(const csr_matrix& x)
      : m_m(x.m_m), m_n(x.m_n), m_nnz(x.m_nnz), m_base(x.m_base),
        m_value(x.m_value), m_rowptr(x.m_rowptr), m_colind(x.m_colind) {}
    virtual ~csr_matrix() {}
    int size() const { return m_nnz; }
    int nrow() const { return m_m; }
    int ncol() const { return m_n; }
    int nnz() const { return m_nnz; }
    int base() const { return m_base; }

    double& operator[](int i) { return m_value[i]; }
    const double& operator[](int i) const { return m_value[i]; }
    const int& rowptr(int i) const { return m_rowptr[i]; }
    const int& colind(int i) const { return m_colind[i]; }
  };

  class csc_matrix {
  private:
    const int m_m;
    const int m_n;
    const int m_nnz;
    const int m_base;

    double* m_value;
    const int* m_colptr;
    const int* m_rowind;

  public:
    csc_matrix(int m, int n, int nnz, int base, double* value, const int* colptr, const int* rowind)
      : m_m(m), m_n(n), m_nnz(nnz), m_base(base),
        m_value(value), m_colptr(colptr), m_rowind(rowind) {}
    csc_matrix(const csc_matrix& x)
      : m_m(x.m_m), m_n(x.m_n), m_nnz(x.m_nnz), m_base(x.m_base),
        m_value(x.m_value), m_colptr(x.m_colptr), m_rowind(x.m_rowind) {}
    virtual ~csc_matrix() {}
    int size() const { return m_nnz; }
    int nrow() const { return m_m; }
    int ncol() const { return m_n; }
    int nnz() const { return m_nnz; }
    int base() const { return m_base; }

    double& operator[](int i) { return m_value[i]; }
    const double& operator[](int i) const { return m_value[i]; }
    const int& colptr(int i) const { return m_colptr[i]; }
    const int& rowind(int i) const { return m_rowind[i]; }
  };

  class coo_matrix {
  private:
    const int m_m;
    const int m_n;
    const int m_nnz;
    const int m_base;

    double* m_value;
    const int* m_rowind;
    const int* m_colind;

  public:
    coo_matrix(int m, int n, int nnz, int base, double* value, const int* rowind, const int* colind)
      : m_m(m), m_n(n), m_nnz(nnz), m_base(base),
        m_value(value), m_rowind(rowind), m_colind(colind) {}
    coo_matrix(const coo_matrix& x)
      : m_m(x.m_m), m_n(x.m_n), m_nnz(x.m_nnz), m_base(x.m_base),
        m_value(x.m_value), m_rowind(x.m_rowind), m_colind(x.m_colind) {}
    virtual ~coo_matrix() {}
    int size() const { return m_nnz; }
    int nrow() const { return m_m; }
    int ncol() const { return m_n; }
    int nnz() const { return m_nnz; }
    int base() const { return m_base; }

    double& operator[](int i) { return m_value[i]; }
    const double& operator[](int i) const { return m_value[i]; }
    const int& rowind(int i) const { return m_rowind[i]; }
    const int& colind(int i) const { return m_colind[i]; }
  };

}
