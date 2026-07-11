#ifndef METRIC_H
#define METRIC_H

#include "math/Matrix.hpp"
#include "math/MatrixGMP.hpp"
#include "math/Vector.hpp"
#include "math/VectorGMP.hpp"

constexpr double M = 1;
constexpr double rs = 2 * M;

Vector<6> getYPrime(Vector<6> y);
Vector<8> getYPrime(Vector<8> y);
Vector<3> getYPrime(Vector<3> y);

Matrix<6, 6> getYPrimePrime(Vector<6> y);
Matrix<3, 3> getYPrimePrime(Vector<3> y);

VectorGMP<6> getYPrime(VectorGMP<6> y);
VectorGMP<8> getYPrime(VectorGMP<8> y);
VectorGMP<3> getYPrime(VectorGMP<3> y);

MatrixGMP<6, 6> getYPrimePrime(VectorGMP<6> y);
MatrixGMP<3, 3> getYPrimePrime(VectorGMP<3> y);

#endif
