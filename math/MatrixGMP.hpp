#ifndef MATRIXGMP_HPP
#define MATRIXGMP_HPP

#include "../mpreal.h"
#include "Vector.hpp"
#include "VectorGMP.hpp"
#include <gmpxx.h>
template <int m, int n> class MatrixGMP {

  mpfr::mpreal val[m][n];

public:
  MatrixGMP() {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        val[i][j] = "0";
      }
    }
  }

  MatrixGMP(mpfr::mpreal (&values)[m][n]) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        val[i][j] = values[i][j];
      }
    }
  }

  mpfr::mpreal operator()(int i, int j) const { return val[i][j]; }

  mpfr::mpreal &operator()(int i, int j) { return val[i][j]; }

  MatrixGMP<n, m> transposed() const {
    MatrixGMP<n, m> out;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        out(j, i) = val[i][j];
      }
    }
    return out;
  }

  const VectorGMP<n> row(int i) const { return VectorGMP<n>(val[i]); }

  void setRow(int i, VectorGMP<n> a) {
    for (int j = 0; j < n; j++) {
      val[i][j] = a(j);
    }
  }

  const VectorGMP<m> col(int i) const {
    return VectorGMP<m>(transposed().row(i));
  }

  template <int mSize, int nSize>
  friend MatrixGMP<mSize, nSize> operator+(const MatrixGMP<mSize, nSize> &a,
                                           const MatrixGMP<mSize, nSize> &b);
  template <int mSize, int nSize>
  friend MatrixGMP<mSize, nSize> operator-(const MatrixGMP<mSize, nSize> &a,
                                           const MatrixGMP<mSize, nSize> &b);
  template <int mSize, int nSize>
  friend MatrixGMP<mSize, nSize> operator*(const mpfr::mpreal &a,
                                           const MatrixGMP<mSize, nSize> &b);
  template <int mSize, int nSize>
  friend MatrixGMP<mSize, nSize> operator*(const MatrixGMP<mSize, nSize> &a,
                                           const mpfr::mpreal &b);
  template <int mSize, int nSize, int lSize>
  friend MatrixGMP<mSize, nSize> operator*(const MatrixGMP<mSize, lSize> &a,
                                           const MatrixGMP<lSize, nSize> &b);
  template <int mSize, int nSize>
  friend VectorGMP<mSize> operator*(const MatrixGMP<mSize, nSize> &a,
                                    const VectorGMP<nSize> &b);

  std::string toString() const {
    std::ostringstream ss;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        ss << " " << this->operator()(i, j);
      }
      ss << "\n";
    }
    return ss.str();
  }
};

template <int n> void LU(MatrixGMP<n, n> &a, Vector<n> &p) {
  for (int j = 0; j < n; j++) {
    p(j) = j;
    mpfr::mpreal alpha = mpfr::abs(a(j, j));
    for (int i = j + 1; i < n; i++) {
      if (mpfr::abs(a(i, j)) > alpha) {
        alpha = mpfr::abs(a(i, j));
        p(j) = i;
      }
    }
    if (p(j) != j) {
      VectorGMP<n> rowJ = a.row(j);
      a.setRow(j, a.row((int)p(j)));
      a.setRow((int)p(j), rowJ);
    }
    for (int i = j + 1; i < n; i++) {
      a(i, j) /= a(j, j);
      for (int l = j + 1; l < n; l++) {
        a(i, l) -= a(i, j) * a(j, l);
      }
    }
  }
}

template <int n> void solve(MatrixGMP<n, n> &a, Vector<n> &p, VectorGMP<n> &b) {
  for (int i = 0; i < n; i++) {
    if (p(i) != i) {
      mpfr::mpreal tmp = b(i);
      b(i) = b((int)p(i));
      b((int)p(i)) = tmp;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      b(i) -= a(i, j) * b(j);
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++) {
      b(i) -= a(i, j) * b(j);
    }
    b(i) /= a(i, i);
  }
}

template <int n> mpfr::mpreal det(MatrixGMP<n, n> a) {
  Vector<n> p;
  LU(a, p);
  mpfr::mpreal det = 1;
  for (int i = 0; i < n; i++) {
    if (p(i) != i) {
      det *= -1;
    }
    det *= a(i, i);
  }
  return det;
}

mpfr::mpreal det(const MatrixGMP<1, 1> &a);

mpfr::mpreal det(const MatrixGMP<2, 2> &a);

mpfr::mpreal det(const MatrixGMP<3, 3> &a);

template <int n> MatrixGMP<n, n> inv(MatrixGMP<n, n> a) {
  if (det(a) == 0) {
    return a;
  }
  Vector<n> p;
  MatrixGMP<n, n> inv;
  LU(a, p);
  for (int i = 0; i < n; i++) {
    VectorGMP<n> b;
    b(i) = 1;
    solve(a, p, b);
    inv.setRow(i, b);
  }
  return inv.transposed();
}

template <int mSize, int nSize>
MatrixGMP<mSize, nSize> operator+(const MatrixGMP<mSize, nSize> &a,
                                  const MatrixGMP<mSize, nSize> &b) {
  MatrixGMP<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a(i, j) + b(i, j);
    }
  }
  return out;
}

template <int mSize, int nSize>
MatrixGMP<mSize, nSize> operator-(const MatrixGMP<mSize, nSize> &a,
                                  const MatrixGMP<mSize, nSize> &b) {
  MatrixGMP<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a(i, j) - b(i, j);
    }
  }
  return out;
}

template <int mSize, int nSize>
MatrixGMP<mSize, nSize> operator*(const MatrixGMP<mSize, nSize> &a,
                                  const mpfr::mpreal &b) {
  MatrixGMP<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a(i, j) * b;
    }
  }
  return out;
}

template <int mSize, int nSize>
MatrixGMP<mSize, nSize> operator*(const mpfr::mpreal &a,
                                  const MatrixGMP<mSize, nSize> &b) {
  MatrixGMP<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a * b(i, j);
    }
  }
  return out;
}

template <int mSize, int nSize, int lSize>
MatrixGMP<mSize, nSize> operator*(const MatrixGMP<mSize, lSize> &a,
                                  const MatrixGMP<lSize, nSize> &b) {
  MatrixGMP<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a.row(i) * b.col(j);
    }
  }
  return out;
}

template <int mSize, int nSize>
VectorGMP<mSize> operator*(const MatrixGMP<mSize, nSize> &a,
                           const VectorGMP<nSize> &b) {
  VectorGMP<mSize> out;
  for (int i = 0; i < mSize; i++) {
    out(i) = b * a.row(i);
  }
  return out;
}
#endif
