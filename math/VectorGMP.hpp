#ifndef VECTORGMP_HPP
#define VECTORGMP_HPP

#include "../mpreal.h"
#include <cmath>
#include <gmpxx.h>
#include <sstream>
template <int size> class VectorGMP {

  mpfr::mpreal val[size];

public:
  VectorGMP() {
    for (int i = 0; i < size; i++) {
      val[i] = "0";
    }
  }

  VectorGMP(const mpfr::mpreal (&values)[size]) {
    for (int i = 0; i < size; i++) {
      val[i] = values[i];
    }
  }

  mpfr::mpreal operator()(int i) const { return val[i]; }

  mpfr::mpreal &operator()(int i) { return val[i]; }

  mpfr::mpreal magnitude() const {
    mpfr::mpreal sum("0");
    for (int i = 0; i < size; i++) {
      sum += (val[i] * val[i]);
    }
    return mpfr::sqrt(sum);
  }

  VectorGMP<size> normalized(mpfr::mpreal a) const {
    VectorGMP<size> norm;
    mpfr::mpreal coefficient = a / magnitude();
    for (int i = 0; i < size; i++) {
      norm(i) = val[i] * coefficient;
    }
    return norm;
  }

  void normalize(mpfr::mpreal a) {
    mpfr::mpreal coefficient = a / magnitude();
    for (int i = 0; i < size; i++) {
      val[i] *= coefficient;
    }
  }

  template <int s>
  friend VectorGMP<s> operator+(const VectorGMP<s> &a, const VectorGMP<s> &b);
  template <int s>
  friend VectorGMP<s> operator-(const VectorGMP<s> &a, const VectorGMP<s> &b);
  template <int s>
  friend VectorGMP<s> operator*(const mpfr::mpreal &a, const VectorGMP<s> &b);
  template <int s>
  friend VectorGMP<s> operator*(const VectorGMP<s> &a, const mpfr::mpreal &b);
  template <int s>
  friend mpfr::mpreal operator*(const VectorGMP<s> &a, const VectorGMP<s> &b);

  std::string toString() {
    std::stringstream stream;
    for (int i = 0; i < size; i++) {
      stream << " " << val[i];
    }
    return stream.str();
  }
};

template <int s>
VectorGMP<s> operator+(const VectorGMP<s> &a, const VectorGMP<s> &b) {
  VectorGMP<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) + b(i);
  }
  return out;
}

template <int s>
VectorGMP<s> operator-(const VectorGMP<s> &a, const VectorGMP<s> &b) {
  VectorGMP<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) - b(i);
  }
  return out;
}

template <int s>
VectorGMP<s> operator*(const mpfr::mpreal &a, const VectorGMP<s> &b) {
  VectorGMP<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a * b(i);
  }
  return out;
}

template <int s>
VectorGMP<s> operator*(const VectorGMP<s> &a, const mpfr::mpreal &b) {
  VectorGMP<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) * b;
  }
  return out;
}

template <int s>
mpfr::mpreal operator*(const VectorGMP<s> &a, const VectorGMP<s> &b) {
  mpfr::mpreal sum = 0;
  for (int i = 0; i < s; i++) {
    sum += a(i) * b(i);
  }
  return sum;
}

#endif
