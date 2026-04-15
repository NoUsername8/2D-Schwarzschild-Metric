#include "GeometryHelper.h"
#include "Metric.h"
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

vec6 getYPrime(vec6 y) {
  double r = y.y1;
  double dr = y.y2;
  double dt = y.y3;
  double h = y.y4;
  double dh = y.y5;
  double dp = y.y6;

  double a = r - rs;
  return vec6{dr,
              -((rs / 2) * a / (r * r * r)) * dt * dt +
                  (rs / 2) / (r * a) * dr * dr + a * dh * dh +
                  a * std::sin(h) * std::sin(h) * dp * dp,
              -rs / (a * r) * dr * dt,
              dh,
              -2 / r * dr * dh + std::sin(h) * std::cos(h) * dp * dp,
              -2 / r * dr * dp - 2 * cot(h) * dh * dp};
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

matrix6x6 getYPrimePrime(vec6 y) {
  double r = y(0);
  double dr = y(1);
  double dt = y(2);
  double h = y(3);
  double dh = y(4);
  double dp = y(5);

  return matrix6x6{
      {0, 1, 0, 0, 0, 0},
      {dp * dp * std::pow(std::sin(h), 2) + dh * dh -
           (rs * (2 * r - 3 * rs)) / (2 * std::pow(r, 4)) * dt * dt +
           (rs * (rs - 2 * r)) / (2 * r * r * std::pow(r - rs, 2)) * dr * dr,
       rs / (r * r - r * rs) * dr, (rs * (r - rs)) / (std::pow(r, 3)) * dt,
       2 * dp * dp * std::sin(h) * std::cos(h) * (r - rs), 2 * dh * (r - rs),
       2 * dp * std::pow(std::sin(h), 2) * (r - rs)},
      {(rs * (2 * r - rs)) / (r * r * std::pow(r - rs, 2)) * dr * dt,
       -rs / (r * (r - rs)) * dt, -rs / (r * (r - rs)) * dr, 0, 0, 0

      },
      {0, 0, 0, 0, 1, 0},
      {2 / (r * r) * dr * dh, -2 / r * dh, 0,
       dp * dp * (std::pow(std::cos(h), 2) - std::pow(std::sin(h), 2)),
       -2 / r * dr, 2 * dp * std::sin(h) * std::cos(h)},
      {2 / (r * r) * dr * dp, -2 / r * dp, 0, 2 * dh * dp * std::pow(csc(h), 2),
       -2 * dp * cot(h), -(2 * (dh * r * cot(h) + dr)) / r}};
}
