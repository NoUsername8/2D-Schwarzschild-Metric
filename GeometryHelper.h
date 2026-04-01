#include <sstream>
#include <string>
#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

//------------------------------------------------------------------------------

struct vec2 {
  union {
    struct {
      double x;
      double y;
    };
    double val[2];
  };

  double operator()(int i) const { return val[i]; }

  double &operator()(int i) { return val[i]; }
};

vec2 operator+(vec2 a, vec2 b);
vec2 operator-(vec2 a, vec2 b);
vec2 operator*(double a, vec2 v);
vec2 operator*(vec2 v, double a);

//------------------------------------------------------------------------------

struct vec6 {
  union {
    struct {
      double y1;
      double y2;
      double y3;
      double y4;
      double y5;
      double y6;
    };
    double y[6];
  };

  double &val(int i) { return y[i]; }

  double operator()(int i) const { return y[i]; }

  double &operator()(int i) { return y[i]; }

  std::string toString() {
    std::stringstream ss;
    for (int i = 0; i < 6; ++i) {
      ss << " " << (int)val(i);
    }
    return ss.str();
  }
};

vec6 operator+(vec6 a, vec6 b);
vec6 operator-(vec6 a, vec6 b);
vec6 operator*(vec6 v, double a);
vec6 operator*(double a, vec6 v);

//------------------------------------------------------------------------------

struct vec8 {
  union {
    struct {
      double y1;
      double y2;
      double y3;
      double y4;
      double y5;
      double y6;
      double y7;
      double y8;
    };
    double y[8];
  };

  double operator()(int i) const { return y[i]; }

  double &operator()(int i) { return y[i]; }
};

vec8 operator+(vec8 a, vec8 b);
vec8 operator-(vec8 a, vec8 b);
vec8 operator*(vec8 v, double a);
vec8 operator*(double a, vec8 v);

//------------------------------------------------------------------------------

struct Basis {
  vec2 e1;
  vec2 e2;
};

struct TransformingVector {
  Basis basis;
  vec2 components;
};

TransformingVector PolarToOrthonormalBasis(TransformingVector v);

double getAngle(double x, double y);

Basis getPolarBasis(double r, double theta);

vec2 CartesianTransformaion(double r, double t);

vec2 PolarTransformation(double x, double y);

TransformingVector UpdateBasis(TransformingVector v, double t, double r);

vec2 CalcPolarPositionFormBasis(Basis b);

double cot(double x);

double csc(double x);

double Dot6(vec6 a, vec6 b);

//------------------------------------------------------------------------------

struct matrix6x6 {
  union {
    struct {
      vec6 row1;
      vec6 row2;
      vec6 row3;
      vec6 row4;
      vec6 row5;
      vec6 row6;
    };
    vec6 row[6];
  };

  matrix6x6() = default;

  matrix6x6(vec6 row1, vec6 row2, vec6 row3, vec6 row4, vec6 row5, vec6 row6)
      : row1{row1}, row2{row2}, row3{row3}, row4{row4}, row5{row5}, row6{row6} {
  }

  vec6 operator*(vec6 a) {
    return vec6{Dot6(row1, a), Dot6(row2, a), Dot6(row3, a),
                Dot6(row4, a), Dot6(row5, a), Dot6(row6, a)};
  }

  vec6 &getRow(int index) { return row[index]; }

  double operator()(int i, int j) const { return row[i](j); }

  double &operator()(int i, int j) { return row[i](j); }

  matrix6x6 transpose() {
    matrix6x6 mat;
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        mat.getRow(i).val(j) = getRow(j).val(i);
      }
    }
    return mat;
  }

  matrix6x6 clone() {
    matrix6x6 mat;
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        mat.getRow(i).val(j) = getRow(i).val(j);
      }
    }
    return mat;
  }

  double determinant() {
    matrix6x6 tmp = clone();
    vec6 p;
    tmp.LU(p);
    double det = 1;
    for (int i = 0; i < 6; i++) {
      if (p.val(i) != i) {
        det *= -1;
      }
      det *= tmp.getRow(i).val(i);
    }
    return det;
  }

  matrix6x6 inverse() {
    if (determinant() == 0) {
      return clone();
    }

    matrix6x6 tmp = clone();
    vec6 p;
    matrix6x6 inv;
    tmp.LU(p);
    for (int i = 0; i < 6; i++) {
      vec6 b{0, 0, 0, 0, 0, 0};
      b.val(i) = 1;

      tmp.solve(p, b);
      inv.getRow(i) = b;
    }
    return inv.transpose();
  }

  void solve(vec6 p, vec6 &b) {
    for (int i = 0; i < 6; i++) {
      if (p.val(i) != i) {
        double tmp = b.val(i);
        b.val(i) = b.val((int)p.val(i));
        b.val((int)p.val(i)) = tmp;
      }
    }
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < i; j++) {
        b.val(i) -= getRow(i).val(j) * b.val(j);
      }
    }
    for (int i = 5; i >= 0; i--) {
      for (int j = i + 1; j < 6; j++) {
        b.val(i) -= getRow(i).val(j) * b.val(j);
      }
      b.val(i) /= getRow(i).val(i);
    }
  }

  void LU(vec6 &p) {

    for (int j = 0; j < 6; j++) {
      p.val(j) = j;
      double alpha = std::abs(getRow(j).val(j));
      for (int i = j + 1; i < 6; i++) {
        if (std::abs(getRow(i).val(j)) > alpha) {
          alpha = std::abs(getRow(i).val(j));
          p.val(j) = i;
        }
      }
      if (p.val(j) != j) {
        vec6 rowJ = getRow(j);
        getRow(j) = getRow((int)p.val(j));
        getRow((int)p.val(j)) = rowJ;
      }
      for (int i = j + 1; i < 6; i++) {
        getRow(i).val(j) /= getRow(j).val(j);
        for (int l = j + 1; l < 6; l++) {
          getRow(i).val(l) -= getRow(i).val(j) * getRow(j).val(l);
        }
      }
    }
  }

  std::string toString() {
    std::ostringstream ss;
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        ss << " " << getRow(i).val(j);
      }
      ss << "\n";
    }
    return ss.str();
  }
};

matrix6x6 id6x6();
matrix6x6 operator*(const matrix6x6 &A, const matrix6x6 &B);
matrix6x6 operator*(double x, const matrix6x6 &A);
matrix6x6 operator*(const matrix6x6 &A, double x);
matrix6x6 operator+(const matrix6x6 &A, const matrix6x6 &B);
matrix6x6 operator-(const matrix6x6 &A, const matrix6x6 &B);

#endif
