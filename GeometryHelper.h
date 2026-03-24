#include <string>
#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

struct vec2 {
  double x;
  double y;

  vec2 operator+(const vec2 a) { return vec2{a.x + x, a.y + y}; }

  vec2 operator-(vec2 a) { return vec2{-a.x + x, -a.y + y}; }
};

struct vec6 {
  double y1;
  double y2;
  double y3;
  double y4;
  double y5;
  double y6;

  vec6 operator+(vec6 a) {
    return vec6{a.y1 + y1, a.y2 + y2, a.y3 + y3,
                a.y4 + y4, a.y5 + y5, a.y6 + y6};
  }

  vec6 operator*(double a) {
    return vec6{a * y1, a * y2, a * y3, a * y4, a * y5, a * y6};
  }

  double* val(int index) {
    switch (index) {
    case 1:
      return &y1;
      break;
    case 2:
      return &y1;
      break;
    case 3:
      return &y1;
      break;
    case 4:
      return &y1;
      break;
    case 5:
      return &y1;
      break;
    case 6:
      return &y1;
      break;
    default:
      return &y1;
      break;
    }
  }
};

struct vec8 {
  double y1;
  double y2;
  double y3;
  double y4;
  double y5;
  double y6;
  double y7;
  double y8;

  vec8 operator+(vec8 a) {
    return vec8{a.y1 + y1, a.y2 + y2, a.y3 + y3, a.y4 + y4,
                a.y5 + y5, a.y6 + y6, a.y7 + y7, a.y8 + y8};
  }

  vec8 operator*(double a) {
    return vec8{a * y1, a * y2, a * y3, a * y4, a * y5, a * y6, a * y7, a * y8};
  }
};

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
double Dot6(vec6 a, vec6 b);

struct matrix6x6{

  vec6 row1;
  vec6 row2;
  vec6 row3;
  vec6 row4;
  vec6 row5;
  vec6 row6;

  vec6 operator*(vec6 a) {
    return vec6{
      Dot6(row1, a),
      Dot6(row2, a),
      Dot6(row3, a),
      Dot6(row4, a),
      Dot6(row5, a),
      Dot6(row6, a)
    };
  }

  vec6* getRow(int index) {
    switch (index) {
      case 1:
        return &row1;
        break;
      case 2:
        return &row2;
        break;
      case 3:
        return &row3;
        break;
      case 4:
        return &row4;
        break;
      case 5:
        return &row5;
        break;
      case 6:
        return &row6;
        break;
      default:
        return &row1;
        break;
    }
  }

  matrix6x6 transpose() {
    matrix6x6 mat;
    for(int i = 0; i < 6; i++) {
      for(int j = 0; j < 6; j++) {
        *(*mat.getRow(i)).val(j) = *(*getRow(j)).val(i);
      } 
    }
    return mat;
  }

  matrix6x6 clone() {
    matrix6x6 mat;
    for(int i = 0; i < 6; i++) {
      for(int j = 0; j < 6; j++) {
        *(*mat.getRow(i)).val(j) = *(*getRow(i)).val(j);
      } 
    }
    return mat;
  }

  double determinant() {
    matrix6x6 tmp = clone();
    vec6 p;
    tmp.LU(p);
	  double det = 1;
    for(int i = 0; i < 6; i++) {
	    if (*p.val(i) != i) {
		    det *= -1;
	    }
	    det *= *(*tmp.getRow(i)).val(i);
    }
    return det;
  }

  matrix6x6 inverse() {
    if(determinant() == 0) {
      return clone();
    }

    matrix6x6 tmp = clone();
    vec6 p;
    matrix6x6 inv;
    tmp.LU(p);
    for(int i = 0; i < 6; i++) {
	    vec6 b{0, 0, 0, 0, 0, 0};
	    *b.val(i) = 1;
	    tmp.solve(p, b);
	    *inv.getRow(i) = b;
	  }
	  return inv.transpose();
  }

  void solve(vec6 p, vec6 b) {
    for(int i = 0; i < 6; i++) {
      if(*p.val(i) != i) {
        double tmp = *b.val(i);
        *b.val(i) = *b.val((int)*p.val(i));
        *b.val((int)*p.val(i)) = tmp;
      }
    }
    for(int i = 0; i < 6; i++) {
      for(int j = 0; j < i; j++) {
        *b.val(i) -= *(*getRow(i)).val(j) * *b.val(j);
      }
    }
    for(int i = 5; i >= 0; i--) {
      for(int j = i + 1; j < 6; j++) {
        *b.val(i) -= *(*getRow(i)).val(j) * *b.val(j);
      }
      *b.val(i) /= *(*getRow(i)).val(i);
    }
  }

  void LU(vec6 p) {;
    for(int j = 0; j < 6; j++) {
      *p.val(j) = j;
      double alpha = std::abs(*(*getRow(j)).val(j));
      for(int i = j + 1; i < 6; i++) {
        if(std::abs(*(*getRow(i)).val(j)) > alpha) {
          alpha = std::abs(*(*getRow(i)).val(j));
          *p.val(j) = i;
        }
      }
      if(*p.val(j) != j) {
        vec6 rowJ = *getRow(j);
        *getRow(j) = *getRow((int)*p.val(j));
        *getRow((int)*p.val(j)) = rowJ;
      }
      for(int i = j + 1; i < 6; i++) {
        *(*getRow(i)).val(j) /= *(*getRow(j)).val(j);
        for(int l = j + 1; l < 6; l++) {
          *(*getRow(i)).val(l) -= *(*getRow(i)).val(j) * *(*getRow(j)).val(l);
        }
      }
    }
  }

  std::string toString() {
    std::string string = "";
    for(int i = 0; i < 6; i++) {
      for(int j = 0; j < 6; j++) {
        string += *(*getRow(i)).val(j);
      } 
      string += "\n";
    }
    return string;
  }
};

#endif
