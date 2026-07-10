#ifndef VECTORGMP_HPP
#define VECTORGMP_HPP

#include <cmath>
#include <gmpxx.h>
#include <sstream>
template <int size> class VectorGMP {

  mpf_class val[size];

public:
  VectorGMP() {
    for (int i = 0; i < size; i++) {
      val[i] = "0";
    }
  }

  VectorGMP(const mpf_class (&values)[size]) {
    for (int i = 0; i < size; i++) {
      val[i] = values[i];
    }
  }

  mpf_class operator()(int i) const { return val[i]; }

  mpf_class &operator()(int i) { return val[i]; }

  mpf_class magnitude() const {
    mpf_class sum("0");
    for (int i = 0; i < size; i++) {
      sum += (val[i] * val[i]);
    }
    return sqrt(sum);
  }

  VectorGMP<size> normalized(mpf_class a) const {
    VectorGMP<size> norm;
    mpf_class coefficient = a / magnitude();
    for (int i = 0; i < size; i++) {
      norm(i) = val[i] * coefficient;
    }
    return norm;
  }

  void normalize(mpf_class a) {
    mpf_class coefficient = a / magnitude();
    for (int i = 0; i < size; i++) {
      val[i] *= coefficient;
    }
  }

  template <int s>
  friend VectorGMP<s> operator+(const VectorGMP<s> &a, const VectorGMP<s> &b);
  template <int s>
  friend VectorGMP<s> operator-(const VectorGMP<s> &a, const VectorGMP<s> &b);
  template <int s>
  friend VectorGMP<s> operator*(const double &a, const VectorGMP<s> &b);
  template <int s>
  friend VectorGMP<s> operator*(const VectorGMP<s> &a, const double &b);
  template <int s>
  friend double operator*(const VectorGMP<s> &a, const VectorGMP<s> &b);

  std::string toString() {
    std::stringstream stream;
    for (int i = 0; i < size; i++) {
      stream << " " << val[i].get_str();
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
VectorGMP<s> operator*(const double &a, const VectorGMP<s> &b) {
  VectorGMP<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a * b(i);
  }
  return out;
}

template <int s>
VectorGMP<s> operator*(const VectorGMP<s> &a, const double &b) {
  VectorGMP<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) * b;
  }
  return out;
}

template <int s>
mpf_class operator*(const VectorGMP<s> &a, const VectorGMP<s> &b) {
  mpf_class sum = 0;
  for (int i = 0; i < s; i++) {
    sum += a(i) * b(i);
  }
  return sum;
}

#endif
