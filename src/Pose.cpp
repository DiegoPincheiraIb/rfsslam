#include "Pose.hpp"

/********** Implementation of example 2d vehicle pose state **********/

Pose2d::Pose2d(){
  StateType x;
  x << 0, 0, 0;
  set(x);
}

Pose2d::Pose2d(StateType x){
  set(x);
}

Pose2d::Pose2d( double x, double y, double theta ){
  StateType state;
  state << x, y, theta;
  set(state);
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

