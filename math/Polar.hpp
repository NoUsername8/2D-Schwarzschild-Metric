#ifndef POLAR_HPP
#define POLAR_HPP

#include "Vector.hpp"
#include "VectorGMP.hpp"

Vector<2> polarToCartesian(Vector<2> polarCoords);
Vector<2> cartesianToPolar(Vector<2> cartesianCoords);
double cot(double x);
VectorGMP<2> polarToCartesian(VectorGMP<2> polarCoords);
VectorGMP<2> cartesianToPolar(VectorGMP<2> cartesianCoords);
#endif
