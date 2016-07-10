#ifndef POINT_H_
#define POINT_H_
#include <assert.h>

const int DIM = 2;


class Point {
   float coord[DIM];
   public:
   Point() {}
   Point(float x, float y) {
      coord[0] = x;
      coord[1] = y;
   }

   Point operator+(Point b);
   Point operator-(Point b);
   Point operator*(float lambda);
   Point operator/(float lambda);
   float operator*(Point b);
   float Norm();
   Point Rotate(float angle);
   float  operator[](int i) const;
   float& operator [](int i );
};



Point Point::operator+(Point b) {
   Point c;
   for (int i=0;i<DIM;i++)
      c.coord[i] = coord[i] + b.coord[i];
   return c;
}

Point Point::operator-(Point b) {
   Point c;
   for (int i=0;i<DIM;i++)
      c.coord[i] = coord[i] - b.coord[i];
   return c;
}

Point Point::operator*(float lambda) {
   Point c;
   for (int i=0;i<DIM;i++)
      c.coord[i] = coord[i] * lambda;
   return c;
}

Point Point::operator/(float lambda) {
   Point c;
   for (int i=0;i<DIM;i++)
      c.coord[i] = coord[i] / lambda;
   return c;
}

float Point::operator*(Point b) {
   float sum = 0;
   for (int i=0;i<DIM;i++)
      sum += coord[i] * b.coord[i];
   return sum;
}

float Point::Norm() {
   float sum = 0;
   for (int i=0;i<DIM;i++)
      sum += coord[i] * coord[i];
   return sqrt(sum);
}

float Point::operator []( int i) const {
   assert ( i>=0 || i < DIM);
   return coord[i]; 
}

float& Point::operator []( int i) {
   assert ( i>=0 || i < DIM);
   return coord[i]; 
}

#endif /* POINT_H_ */
