#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <sstream>
template <int size> class Vector {

  double val[size];

public:
  Vector() {
    for (int i = 0; i < size; i++) {
      val[i] = 0;
    }
  }

  Vector(const double (&values)[size]) {
    for (int i = 0; i < size; i++) {
      val[i] = values[i];
    }
  }

  double operator()(int i) const { return val[i]; }

  double &operator()(int i) { return val[i]; }

  double magnitude() const {
    double sum = 0;
    for (int i = 0; i < size; i++) {
      sum += (val[i] * val[i]);
    }
    return std::sqrt(sum);
  }

  Vector<size> normalized(double a) const {
    Vector<size> norm;
    double coefficient = a / magnitude();
    for (int i = 0; i < size; i++) {
      norm(i) = val[i] * coefficient;
    }
    return norm;
  }

  void normalize(double a) {
    double coefficient = a / magnitude();
    for (int i = 0; i < size; i++) {
      val[i] *= coefficient;
    }
  }

  template <int s>
  friend Vector<s> operator+(const Vector<s> &a, const Vector<s> &b);
  template <int s>
  friend Vector<s> operator-(const Vector<s> &a, const Vector<s> &b);
  template <int s>
  friend Vector<s> operator*(const double &a, const Vector<s> &b);
  template <int s>
  friend Vector<s> operator*(const Vector<s> &a, const double &b);
  template <int s>
  friend double operator*(const Vector<s> &a, const Vector<s> &b);

  std::string toString() {
    std::stringstream stream;
    for (int i = 0; i < size; i++) {
      stream << " " << val[i];
    }
    return stream.str();
  }
};

template <int s> Vector<s> operator+(const Vector<s> &a, const Vector<s> &b) {
  Vector<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) + b(i);
  }
  return out;
}

template <int s> Vector<s> operator-(const Vector<s> &a, const Vector<s> &b) {
  Vector<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) - b(i);
  }
  return out;
}

template <int s> Vector<s> operator*(const double &a, const Vector<s> &b) {
  Vector<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a * b(i);
  }
  return out;
}

template <int s> Vector<s> operator*(const Vector<s> &a, const double &b) {
  Vector<s> out;
  for (int i = 0; i < s; i++) {
    out(i) = a(i) * b;
  }
  return out;
}

template <int s> double operator*(const Vector<s> &a, const Vector<s> &b) {
  double sum = 0;
  for (int i = 0; i < s; i++) {
    sum += a(i) * b(i);
  }
  return sum;
}

#endif
