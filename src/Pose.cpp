#include "Pose.hpp"

/********** Implementation of example 2d vehicle pose state **********/

Pose2d::Pose2d(){}

Pose2d::Pose2d(Vec &x, double t) : 
  PoseWithUncertainty<Vec, Mat>(x, (Mat() << Mat::Zero()).finished() ,t) {}

Pose2d::Pose2d(Vec &x, Mat &Sx, double t) : 
  PoseWithUncertainty<Vec, Mat>(x, Sx, t) {}

Pose2d::Pose2d( double x, double y, double theta, 
		double var_x, double var_y, double var_theta,
		double t ){
  Vec state;
  state << x, y, theta;
  Mat cov;
  cov << 
    var_x, 0, 0,
    0, var_y, 0,
    0, 0, var_x;
  set(state, cov, t);
}

Pose2d::~Pose2d(){}

