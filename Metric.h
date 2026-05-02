#ifndef METRIC_H
#define METRIC_H

#include "GeometryHelper.h"

constexpr double M = 2;
constexpr double rs = 2 * M;

vec6 getYPrime(vec6 y);
vec8 getYPrime(vec8 y);
vec3 getYPrime(vec3 y);

matrix6x6 getYPrimePrime(vec6 y);
matrix3x3 getYPrimePrime(vec3 y);

#endif
