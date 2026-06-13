#include "ScreenHelper.h"
#include "math/Vector.hpp"
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

#define SCREEN_WIDTH 1280
#define SCREEN_HEIGHT 720
#define SCALE 50

Vector<2> TransformToScreenCoords(double x, double y, Vector<2> windowPos) {
  Vector<2> vector;
  vector(0) = (SCREEN_WIDTH / 2.0) - (x * SCALE) + windowPos(0);
  vector(1) = (SCREEN_HEIGHT / 2.0) - (y * SCALE) + windowPos(1);
  return vector;
}

Vector<2> TransformToSimulationCoords(double x, double y, Vector<2> windowPos) {
  Vector<2> vector;
  vector(0) = ((SCREEN_WIDTH / 2.0) - (x - windowPos(0))) / SCALE;
  vector(1) = ((SCREEN_HEIGHT / 2.0) - (y - windowPos(1))) / SCALE;
  return vector;
}
