#include "Matrix.hpp"

double det(const Matrix<1, 1> &a) { return a(0, 0); }

double det(const Matrix<2, 2> &a) {
  return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
}

double det(const Matrix<3, 3> &a) {
  return a(0, 0) * a(1, 1) * a(2, 2) + a(0, 1) * a(1, 2) * a(2, 0) +
         a(0, 2) * a(1, 0) * a(2, 1) - a(0, 0) * a(1, 2) * a(2, 1) -
         a(0, 1) * a(1, 0) * a(2, 2) - a(0, 2) * a(1, 1) * a(1, 0);
}
