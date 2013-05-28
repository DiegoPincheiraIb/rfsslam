#include "Landmark.hpp"

Landmark2d::Landmark2d(){
  x_ << 0, 0;
}

Landmark2d::Landmark2d(vec x){
  setState(x);
}

Landmark2d::~Landmark2d(){}

void Landmark2d::setState(vec x){
  x_ = x;
}

void Landmark2d::getState(vec& x){
  x = x_;
}



Landmark1d::Landmark1d(){
  x_ = 0;
}

Landmark1d::Landmark1d(double x){
  setState(x);
}

Landmark1d::~Landmark1d(){}

void Landmark1d::setState(double x){
  x_ = x;
}

void Landmark1d::getState(double& x){
  x = x_;
}
