#include "Pose.hpp"

/********** Implementation of example 2d vehicle pose state **********/

Pose2d::Pose2d(){}

Pose2d::Pose2d(StateType x, double t) : Pose<StateType>(x, t) {}

Pose2d::Pose2d( double x, double y, double theta, double t ){
  StateType state;
  state << x, y, theta;
  set(state, t);
}

Pose2d::~Pose2d(){}

