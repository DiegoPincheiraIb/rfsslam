#include "Pose.hpp"

/********** Implementation of example 2d vehicle pose state **********/

Pose2d::Pose2d(){
  vec x;
  x << 0, 0, 0;
  set(x);
}

Pose2d::Pose2d(vec x){
  set(x);
}

Pose2d::~Pose2d(){}


/********** Implementation of example 1d vehicle pose state **********/

Pose1d::Pose1d(){
  set(0);
}

Pose1d::Pose1d(double x){
  set(x);
}

Pose1d::~Pose1d(){}

