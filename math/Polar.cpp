#include "Polar.hpp"
#include "Vector.hpp"
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
