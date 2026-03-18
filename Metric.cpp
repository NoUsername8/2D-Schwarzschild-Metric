#include "GeometryHelper.h"
#include "Metric.h"
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

vec6 getYPrime(vec6 y) {
  double a = y.y1 - rs;
  return vec6{y.y2,
              -((rs / 2) * a / (y.y1 * y.y1 * y.y1)) * y.y3 * y.y3 +
                  (rs / 2) / (y.y1 * a) * y.y2 * y.y2 + a * y.y5 * y.y5 +
                  a * std::sin(y.y4) * std::sin(y.y4) * y.y6 * y.y6,
              -rs / (a * y.y1) * y.y2 * y.y3,
              y.y5,
              -2 / y.y1 * y.y2 * y.y5 +
                  std::sin(y.y4) * std::cos(y.y4) * y.y6 * y.y6,
              -2 / y.y1 * y.y2 * y.y6 - 2 * cot(y.y4) * y.y5 * y.y6};
}

vec8 getYPrime(vec8 y) {
  double a = y.y1 - rs;
  return vec8{y.y2,
              -((rs / 2) * a / (y.y1 * y.y1 * y.y1)) * y.y4 * y.y4 +
                  (rs / 2) / (y.y1 * a) * y.y2 * y.y2 + a * y.y6 * y.y6 +
                  a * std::sin(y.y5) * std::sin(y.y5) * y.y8 * y.y8,
              y.y4,
              -rs / (a * y.y1) * y.y2 * y.y4,
              y.y6,
              -2 / y.y1 * y.y2 * y.y6 +
                  std::sin(y.y5) * std::cos(y.y5) * y.y8 * y.y8,
              y.y8,
              -2 / y.y1 * y.y2 * y.y8 - 2 * cot(y.y5) * y.y6 * y.y8};
}
