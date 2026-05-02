#include "GeometryHelper.h"
#include "Metric.h"
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

vec6 getYPrime(vec6 y) {
  double r = y.y1;
  double vr = y.y2;
  double t = y.y3;
  double vt = y.y4;
  double phi = y.y5;
  double vphi = y.y6;

  return vec6{vr,
              -M * (r - 2 * M) * vt * vt / (r * r * r) +
                  M * vr * vr / (r * (r - 2 * M)) + (r - 2 * M) * vphi * vphi,
              vt,
              -2 * M * vr * vt / (r * (r - 2 * M)),
              vphi,
              -2 * vr * vphi / r};
}
