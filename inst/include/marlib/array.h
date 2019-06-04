#pragma once

namespace marlib {

class array {
private:
  int m_n;
  double* m_ptr;
  int m_inc;

public:
  array(int n, double* ptr, int inc) : m_n(n), m_ptr(ptr), m_inc(inc) {}
  virtual ~array() {}
  double& operator[](int i) { return m_ptr[i]; }
  const double& operator[](int i) const { return m_ptr[i]; }
  int size() const { return m_n; }
  int inc() const { return m_inc; }
};

// class array1 {
// private:
//   int m_n;
//   double* m_ptr;
// public:
//   array1(int n) : m_n(n) {
//     m_ptr = new double [m_n];
//   }
//   virtual ~array1() {
//     delete [] m_ptr;
//   }
//   double& operator[](int i) { return m_ptr[i]; }
//   const double& operator[](int i) const { return m_ptr[i]; }
//   int size() const { return m_n; }
// };
//
// class array2 {
// private:
//   int m_m;
//   int m_n;
//   double* m_ptr;
//   std::vector<array*> m_elem;
//
// public:
//   array2(int m, int n) : m_m(m), m_n(n), m_elem(m) {
//     m_ptr = new double [m_m * m_n];
//     double* p = m_ptr;
//     for (int i=0; i<m_m; i++, p+=m_n) {
//       m_elem[i] = new array (m_n, p);
//     }
//   }
//
//   virtual ~array2() {
//     for (int i=0; i<m_m; i++) {
//       delete m_elem[i];
//     }
//     delete [] m_ptr;
//   }
//
//   array& operator[](int i) {
//     return *m_elem[i];
//   }
//
//   const array& operator[](int i) const {
//     return *m_elem[i];
//   }
// };

}
