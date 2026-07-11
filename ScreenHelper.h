#ifndef SCREENHELPER_H
#define SCREENHELPER_H
#include "math/Vector.hpp"
#include "math/VectorGMP.hpp"
#include "mpreal.h"

Vector<2> TransformToScreenCoords(double x, double y, Vector<2> windowPos);
Vector<2> TransformToSimulationCoords(double x, double y, Vector<2> windowPos);

VectorGMP<2> TransformToScreenCoords(mpfr::mpreal x, mpfr::mpreal y,
                                     VectorGMP<2> windowPos);
VectorGMP<2> TransformToSimulationCoords(mpfr::mpreal x, mpfr::mpreal y,
                                         VectorGMP<2> windowPos);

#endif
