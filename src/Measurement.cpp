#include "Measurement.hpp"

/********** Implementation of example 1d measurement **********/

Measurement1d::Measurement1d()
  : Measurement< Eigen::Matrix<double, 1, 1>,
		 Eigen::Matrix<double, 1, 1> > (1)
{}

Measurement1d::Measurement1d( Eigen::Matrix<double, 1, 1> z,
			      Eigen::Matrix<double, 1, 1> Sz, 
			      double t) 
  : Measurement< Eigen::Matrix<double, 1, 1>,
		 Eigen::Matrix<double, 1, 1> > (1, z, Sz, t)
{}

Measurement1d::~Measurement1d(){}

/********** Implementation of example 2d measurement **********/

Measurement2d::Measurement2d() : Measurement<Eigen::Vector2d, Eigen::Matrix2d> (2)
{}

Measurement2d::~Measurement2d(){}

/********** Implementation of example 2d odometry measurement **********/

Odometry2d::Odometry2d() : Measurement< Eigen::Vector3d, Eigen::Matrix3d > (3){}

Odometry2d::Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km, 
		       double vardx_k_km, double vardy_k_km, 
		       double vartheta_k_km, double t) 
  : Measurement< Eigen::Vector3d, Eigen::Matrix3d > (3){
  Eigen::Vector3d u;
  Eigen::Vector3d covu_diag; 
  Eigen::Matrix3d covu;
  u << dx_k_km, dy_k_km, dtheta_k_km;
  covu_diag << vardx_k_km, vardy_k_km, vartheta_k_km;
  covu = covu_diag.asDiagonal();

  set( u, covu, t);
}

Odometry2d::~Odometry2d(){}









