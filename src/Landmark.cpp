#include "Landmark.hpp"

/********** Implementation of example 1d vehicle pose state **********/

Landmark1d::Landmark1d(){
  set(0, 0);
}

Landmark1d::~Landmark1d(){}

double Landmark1d::mahalanobisDist2(double &x){


  double e = x_ - x;

  if(this->pSxInv_ == NULL){
    this->SxInv_ = 1 / this->Sx_;
  }

  return (e * SxInv_ * e);
}

/********** Implementation of example 2d vehicle pose state **********/

Landmark2d::Landmark2d(){
  Eigen::Vector2d x = Eigen::Vector2d::Zero();
  Eigen::Matrix2d Sx = Eigen::Matrix2d::Zero();
  set(x,Sx);
}

Landmark2d::~Landmark2d(){}

double Landmark2d::mahalanobisDist2(Eigen::Vector2d &x){

  Eigen::Vector2d e = x_ - x;

  if(this->pSxInv_ == NULL){
    this->SxInv_ = this->Sx_.inverse();
  }

  return (e.transpose() * SxInv_ * e);
}




