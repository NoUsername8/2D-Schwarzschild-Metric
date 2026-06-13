#ifndef SCREENHELPER_H
#define SCREENHELPER_H
#include "math/Vector.hpp"

Vector<2> TransformToScreenCoords(double x, double y, Vector<2> windowPos);
Vector<2> TransformToSimulationCoords(double x, double y, Vector<2> windowPos);

#endif
