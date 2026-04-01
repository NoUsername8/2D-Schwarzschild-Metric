#ifndef METRIC_H
#define METRIC_H

#include "GeometryHelper.h"

constexpr double M = 2;
constexpr double c = 2;
constexpr double G = 1;
constexpr double rs = (2 * G * M) / (c * c);

vec6 getYPrime(vec6 y);
vec8 getYPrime(vec8 y);

matrix6x6 getYPrimePrime(vec6 y);

#endif
