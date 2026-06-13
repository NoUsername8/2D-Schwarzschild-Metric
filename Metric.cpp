#include "Metric.h"
#include "math/Matrix.hpp"
#include "math/Vector.hpp"
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

Vector<6> getYPrime(Vector<6> y) {
  double r = y(0);
  double vr = y(1);
  double t = y(2);
  double vt = y(3);
  double phi = y(4);
  double vphi = y(5);

  double vals[6] = {vr,
                    -M * (r - 2 * M) * vt * vt / (r * r * r) +
                        M * vr * vr / (r * (r - 2 * M)) +
                        (r - 2 * M) * vphi * vphi,
                    vt,
                    -2 * M * vr * vt / (r * (r - 2 * M)),
                    vphi,
                    -2 * vr * vphi / r};
  return Vector<6>(vals);
}

Vector<3> getYPrime(Vector<3> y) {
  double r = y(0);
  double dr = y(1);
  double dp = y(2);

  double vals[3] = {dr, M * dr * dr / (r * (r - 2 * M)) + (r - 2 * M) * dp * dp,
                    -2 * dr * dp / r};
  return Vector<3>(vals);
}

Matrix<3, 3> getYPrimePrime(Vector<3> y) {
  double r = y(0);
  double dr = y(1);
  double dp = y(2);

  double vals[3][3] = {
      {0, 1, 0},
      {-2 * M * (r - M) * dr * dr / std::pow(r * r - 2 * M * r, 2) + dp * dp,
       2 * M * dr / (r * (r - 2 * M)), 2 * (r - 2 * M) * dp},
      {2 * dr * dp / (r * r), -2 * dp / r, -2 * dr / r}};

  return Matrix<3, 3>(vals);
}
