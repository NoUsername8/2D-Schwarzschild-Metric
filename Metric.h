#ifndef METRIC_H
#define METRIC_H

#include "math/Matrix.hpp"
#include "math/Vector.hpp"

constexpr double M = 1;
constexpr double rs = 2 * M;

Vector<6> getYPrime(Vector<6> y);
Vector<8> getYPrime(Vector<8> y);
Vector<3> getYPrime(Vector<3> y);

Matrix<6, 6> getYPrimePrime(Vector<6> y);
Matrix<3, 3> getYPrimePrime(Vector<3> y);

#endif
