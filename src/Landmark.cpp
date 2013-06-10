#include "Landmark.hpp"

/********** Implementation of example 1d vehicle pose state **********/

Landmark1d::Landmark1d(){
  set(0,0);
}

Landmark1d::~Landmark1d(){}

double Landmark1d::mahalanobisDist(double &x){
  double e = x_ - x;
  return e / Sx_ * e;
}

/********** Implementation of example 2d vehicle pose state **********/

Landmark2d::Landmark2d(){
  Eigen::Vector2d x;
  x << 0, 0;
  Eigen::Matrix2d Sx;
  Sx << 1,0,
        0,1;
  set(x,Sx);
}

Landmark2d::~Landmark2d(){}

double Landmark2d::mahalanobisDist(Eigen::Vector2d &x){
  Eigen::Vector2d e = x_ - x;
  return e.transpose() * Sx_.inverse() * e;
}




