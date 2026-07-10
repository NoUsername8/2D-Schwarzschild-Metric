#include "MatrixGMP.hpp"
#include <gmpxx.h>

mpf_class det(const MatrixGMP<1, 1> &a) { return a(0, 0); }

mpf_class det(const MatrixGMP<2, 2> &a) {
  return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
}

mpf_class det(const MatrixGMP<3, 3> &a) {
  return a(0, 0) * a(1, 1) * a(2, 2) + a(0, 1) * a(1, 2) * a(2, 0) +
         a(0, 2) * a(1, 0) * a(2, 1) - a(0, 0) * a(1, 2) * a(2, 1) -
         a(0, 1) * a(1, 0) * a(2, 2) - a(0, 2) * a(1, 1) * a(1, 0);
}
