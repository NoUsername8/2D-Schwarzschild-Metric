#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"
template <int m, int n> class Matrix {

  double val[m][n];

public:
  Matrix() {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        val[i][j] = 0;
      }
    }
  }

  Matrix(double (&values)[m][n]) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        val[i][j] = values[i][j];
      }
    }
  }

  double operator()(int i, int j) const { return val[i][j]; }

  double &operator()(int i, int j) { return val[i][j]; }

  Matrix<n, m> transposed() const {
    Matrix<n, m> out;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        out(j, i) = val[i][j];
      }
    }
    return out;
  }

  const Vector<n> row(int i) const { return Vector<n>(val[i]); }

  void setRow(int i, Vector<n> a) {
    for (int j = 0; j < n; j++) {
      val[i][j] = a(j);
    }
  }

  const Vector<m> col(int i) const { return Vector<m>(transposed().row(i)); }

  template <int mSize, int nSize>
  friend Matrix<mSize, nSize> operator+(const Matrix<mSize, nSize> &a,
                                        const Matrix<mSize, nSize> &b);
  template <int mSize, int nSize>
  friend Matrix<mSize, nSize> operator-(const Matrix<mSize, nSize> &a,
                                        const Matrix<mSize, nSize> &b);
  template <int mSize, int nSize>
  friend Matrix<mSize, nSize> operator*(const double &a,
                                        const Matrix<mSize, nSize> &b);
  template <int mSize, int nSize>
  friend Matrix<mSize, nSize> operator*(const Matrix<mSize, nSize> &a,
                                        const double &b);
  template <int mSize, int nSize, int lSize>
  friend Matrix<mSize, nSize> operator*(const Matrix<mSize, lSize> &a,
                                        const Matrix<lSize, nSize> &b);
  template <int mSize, int nSize>
  friend Vector<mSize> operator*(const Matrix<mSize, nSize> &a,
                                 const Vector<nSize> &b);

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

template <int n> void LU(Matrix<n, n> &a, Vector<n> &p) {
  for (int j = 0; j < n; j++) {
    p(j) = j;
    double alpha = std::abs(a(j, j));
    for (int i = j + 1; i < n; i++) {
      if (std::abs(a(i, j)) > alpha) {
        alpha = std::abs(a(i, j));
        p(j) = i;
      }
    }
    if (p(j) != j) {
      Vector<n> rowJ = a.row(j);
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

template <int n> void solve(Matrix<n, n> &a, Vector<n> &p, Vector<n> &b) {
  for (int i = 0; i < n; i++) {
    if (p(i) != i) {
      double tmp = b(i);
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

template <int n> double det(Matrix<n, n> a) {
  Vector<n> p;
  LU(a, p);
  double det = 1;
  for (int i = 0; i < n; i++) {
    if (p(i) != i) {
      det *= -1;
    }
    det *= a(i, i);
  }
  return det;
}

double det(const Matrix<1, 1> &a);

double det(const Matrix<2, 2> &a);

double det(const Matrix<3, 3> &a);

template <int n> Matrix<n, n> inv(Matrix<n, n> a) {
  if (det(a) == 0) {
    return a;
  }
  Vector<n> p;
  Matrix<n, n> inv;
  LU(a, p);
  for (int i = 0; i < n; i++) {
    Vector<n> b;
    b(i) = 1;
    solve(a, p, b);
    inv.setRow(i, b);
  }
  return inv.transposed();
}

template <int mSize, int nSize>
Matrix<mSize, nSize> operator+(const Matrix<mSize, nSize> &a,
                               const Matrix<mSize, nSize> &b) {
  Matrix<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a(i, j) + b(i, j);
    }
  }
  return out;
}

template <int mSize, int nSize>
Matrix<mSize, nSize> operator-(const Matrix<mSize, nSize> &a,
                               const Matrix<mSize, nSize> &b) {
  Matrix<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a(i, j) - b(i, j);
    }
  }
  return out;
}

template <int mSize, int nSize>
Matrix<mSize, nSize> operator*(const Matrix<mSize, nSize> &a, const double &b) {
  Matrix<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a(i, j) * b;
    }
  }
  return out;
}

template <int mSize, int nSize>
Matrix<mSize, nSize> operator*(const double &a, const Matrix<mSize, nSize> &b) {
  Matrix<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a * b(i, j);
    }
  }
  return out;
}

template <int mSize, int nSize, int lSize>
Matrix<mSize, nSize> operator*(const Matrix<mSize, lSize> &a,
                               const Matrix<lSize, nSize> &b) {
  Matrix<mSize, nSize> out;
  for (int i = 0; i < mSize; i++) {
    for (int j = 0; j < nSize; j++) {
      out(i, j) = a.row(i) * b.col(j);
    }
  }
  return out;
}

template <int mSize, int nSize>
Vector<mSize> operator*(const Matrix<mSize, nSize> &a, const Vector<nSize> &b) {
  Vector<mSize> out;
  for (int i = 0; i < mSize; i++) {
    out(i) = b * a.row(i);
  }
  return out;
}
#endif
