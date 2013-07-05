#include "Measurement.hpp"

/********** Implementation of example 2d odometry measurement **********/

Odometry2d::Odometry2d(){}

Odometry2d::Odometry2d(Vec &x, Mat &Sx, double t) : 
  RandomVec< Eigen::Vector3d, Eigen::Matrix3d >(x, Sx, t){}

Odometry2d::Odometry2d(Vec &x, double t) : 
  RandomVec< Eigen::Vector3d, Eigen::Matrix3d >(x, t){}

Odometry2d::Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km, 
		       double vardx_k_km, double vardy_k_km, 
		       double vartheta_k_km, double t) {
  Eigen::Vector3d u;
  Eigen::Vector3d covu_diag; 
  Eigen::Matrix3d covu;
  u << dx_k_km, dy_k_km, dtheta_k_km;
  covu_diag << vardx_k_km, vardy_k_km, vartheta_k_km;
  covu = covu_diag.asDiagonal();

  set( u, covu, t);
}

Odometry2d::~Odometry2d(){}









