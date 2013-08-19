#include "Pose.hpp"

/********** Implementation of 1d vechile position state ***********/

Pose1d::Pose1d(){}

Pose1d::Pose1d(double x, double Sx, double t){
  Vec state;
  Mat cov;
  state << x;
  cov << Sx;
  set(state, cov, t);
}

Pose1d::Pose1d(Eigen::Matrix<double, 1, 1> &x, Eigen::Matrix<double, 1, 1> &Sx, double t):
  RandomVec< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >(x, Sx, t){}

Pose1d::Pose1d(double x, double t){
  Vec state;
  state << x;
  set(state, t);
}

Pose1d::Pose1d(Eigen::Matrix<double, 1, 1> &x, double t):
  RandomVec< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >(x, t){}


/********** Implementation of 2d vehicle pose state **********/

Pose2d::Pose2d(){}

Pose2d::Pose2d(Vec &x, Mat &Sx, double t) :
  RandomVec< Eigen::Vector3d, Eigen::Matrix3d >(x, Sx, t){}

Pose2d::Pose2d(Vec &x, double t) : 
  RandomVec< Eigen::Vector3d, Eigen::Matrix3d >(x, t){}

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

