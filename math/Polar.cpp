#include "Polar.hpp"
#include "Vector.hpp"
#include "VectorGMP.hpp"
#include <cmath>

Vector<2> polarToCartesian(Vector<2> polarCoords) {
  double val[2] = {polarCoords(0) * std::cos(polarCoords(1)),
                   polarCoords(0) * std::sin(polarCoords(1))};
  return Vector<2>(val);
}

Vector<2> cartesianToPolar(Vector<2> cartesianCoords) {
  double val[2] = {std::sqrt(std::pow(cartesianCoords(0), 2) +
                             std::pow(cartesianCoords(1), 2)),
                   std::atan2(cartesianCoords(1), cartesianCoords(0))};
  return Vector<2>(val);
}

double cot(double x) { return std::cos(x) == 0 ? 1000000 : 1 / std::cos(x); }

VectorGMP<2> polarToCartesian(VectorGMP<2> polarCoords) {
  mpfr::mpreal val[2] = {polarCoords(0) * mpfr::cos(polarCoords(1)),
                         polarCoords(0) * mpfr::sin(polarCoords(1))};
  return VectorGMP<2>(val);
}

VectorGMP<2> cartesianToPolar(VectorGMP<2> cartesianCoords) {
  mpfr::mpreal val[2] = {mpfr::sqrt(mpfr::pow(cartesianCoords(0), 2) +
                                    mpfr::pow(cartesianCoords(1), 2)),
                         mpfr::atan2(cartesianCoords(1), cartesianCoords(0))};
  return VectorGMP<2>(val);
}
