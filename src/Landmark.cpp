#include "Landmark.hpp"

/********** Implementation of example 1d vehicle pose state **********/

Landmark1d::Landmark1d() 
  : Landmark< Eigen::Matrix<double, 1, 1>,
	      Eigen::Matrix<double, 1, 1> >(1){
  Eigen::Matrix<double, 1, 1> x;
  Eigen::Matrix<double, 1, 1> Sx;
  x << 0;
  Sx << 0;
  set(x, Sx);
}

Landmark1d::~Landmark1d(){}

/********** Implementation of example 2d vehicle pose state **********/

Landmark2d::Landmark2d() : Landmark<Eigen::Vector2d, Eigen::Matrix2d>(3){
  Eigen::Vector2d x = Eigen::Vector2d::Zero();
  Eigen::Matrix2d Sx = Eigen::Matrix2d::Zero();
  set(x,Sx);
}

Landmark2d::~Landmark2d(){}




