#include "Metric.h"
#include "math/Matrix.hpp"
#include "math/Polar.hpp"
#include "math/Vector.hpp"
#include "mpreal.h"
#include <gmpxx.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

//=================================== 3D ======================================
Vector<8> getYPrime(Vector<8> y) {
  double r = y(0);
  double dr = y(1);
  double t = y(2);
  double dt = y(3);
  double h = y(4);
  double dh = y(5);
  double phi = y(6);
  double dphi = y(7);

  double vals[8] = {dr,
                    -M * (r - 2 * M) * dt * dt / (r * r * r) +
                        M * dr * dr / (r * r - 2 * M * r) +
                        (r - 2 * M) * dh * dh +
                        (r - 2 * M) * std::pow(std::sin(h), 2) * dphi * dphi,
                    dt,
                    -2 * M * dr * dt / (r * r - 2 * M * r),
                    dh,
                    -2 * dr * dh / r + std::sin(h) * std::cos(h) * dphi * dphi,
                    dphi,
                    -2 * dr * dh / r - 2 * cot(h) * dh * dphi};
  return Vector<8>(vals);
}
//=================================== 3D mit theta = pi/2 =====================
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

//=================================== 2D ======================================
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

  double vals[3][3] = {{0, 1, 0},
                       {dp * dp - 2 * dr * dr * M * (r - M) /
                                      (r * r * (r - 2 * M) * (r - 2 * M)),
                        2 * M * dr / (r * r - 2 * r * M), 2 * dp * (r - 2 * M)},
                       {2 * dr * dp / (r * r), -2 * dp / r, -2 * dr / r}};

  return Matrix<3, 3>(vals);
}

//=================================== GMP =====================================
//=================================== 3D ======================================
VectorGMP<8> getYPrime(VectorGMP<8> y) {
  mpfr::mpreal r = y(0);
  mpfr::mpreal dr = y(1);
  mpfr::mpreal t = y(2);
  mpfr::mpreal dt = y(3);
  mpfr::mpreal h = y(4);
  mpfr::mpreal dh = y(5);
  mpfr::mpreal phi = y(6);
  mpfr::mpreal dphi = y(7);

  mpfr::mpreal vals[8] = {
      dr,
      -M * (r - 2 * M) * dt * dt / (r * r * r) +
          M * dr * dr / (r * r - 2 * M * r) + (r - 2 * M) * dh * dh +
          (r - 2 * M) * mpfr::pow(mpfr::sin(h), 2) * dphi * dphi,
      dt,
      -2 * M * dr * dt / (r * r - 2 * M * r),
      dh,
      -2 * dr * dh / r + mpfr::sin(h) * mpfr::cos(h) * dphi * dphi,
      dphi,
      -2 * dr * dh / r - 2 * mpfr::cot(h) * dh * dphi};
  return VectorGMP<8>(vals);
}
//=================================== 3D mit theta = pi/2 =====================
VectorGMP<6> getYPrime(VectorGMP<6> y) {
  mpfr::mpreal r = y(0);
  mpfr::mpreal vr = y(1);
  mpfr::mpreal t = y(2);
  mpfr::mpreal vt = y(3);
  mpfr::mpreal phi = y(4);
  mpfr::mpreal vphi = y(5);

  mpfr::mpreal vals[6] = {vr,
                          -M * (r - 2 * M) * vt * vt / (r * r * r) +
                              M * vr * vr / (r * (r - 2 * M)) +
                              (r - 2 * M) * vphi * vphi,
                          vt,
                          -2 * M * vr * vt / (r * (r - 2 * M)),
                          vphi,
                          -2 * vr * vphi / r};
  return VectorGMP<6>(vals);
}

//=================================== 2D ======================================
VectorGMP<3> getYPrime(VectorGMP<3> y) {
  mpfr::mpreal r = y(0);
  mpfr::mpreal dr = y(1);
  mpfr::mpreal dp = y(2);

  mpfr::mpreal vals[3] = {
      dr, M * dr * dr / (r * (r - 2 * M)) + (r - 2 * M) * dp * dp,
      -2 * dr * dp / r};
  return VectorGMP<3>(vals);
}

MatrixGMP<3, 3> getYPrimePrime(VectorGMP<3> y) {
  mpfr::mpreal r = y(0);
  mpfr::mpreal dr = y(1);
  mpfr::mpreal dp = y(2);

  mpfr::mpreal vals[3][3] = {{0, 1, 0},
                             {dp * dp - 2 * dr * dr * M * (r - M) /
                                            (r * r * (r - 2 * M) * (r - 2 * M)),
                              2 * M * dr / (r * r - 2 * r * M),
                              2 * dp * (r - 2 * M)},
                             {2 * dr * dp / (r * r), -2 * dp / r, -2 * dr / r}};

  return MatrixGMP<3, 3>(vals);
}
