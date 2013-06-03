#include <math.h>
#include "ProcessModel.hpp"

/********** Implementation for an examples 2D motion model **********/

Odometry2d::Odometry2d(){}

Odometry2d::Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km, 
		       double vardx_k_km, double vardy_k_km, 
		       double vartheta_k_km, double t){
  Eigen::Vector3d u;
  Eigen::Vector3d covu_diag; 
  Eigen::Matrix3d covu;
  u << dx_k_km, dy_k_km, dtheta_k_km;
  covu_diag << vardx_k_km, vardy_k_km, vartheta_k_km;
  covu = covu_diag.asDiagonal();

  set( u, covu, t);
}

Odometry2d::~Odometry2d(){}



OdometryMotionModel2d::OdometryMotionModel2d(){}

OdometryMotionModel2d::~OdometryMotionModel2d(){}

void OdometryMotionModel2d::step(  Pose2d &s_k, 
				   Pose2d const &s_km, 
				   Odometry2d const &input_k, 
				   double const dT ){
  
  /* State at k-1 
  s_km.get(x_km_i_);
  p_km_i_ = x_km_i_.head(2);
  theta_km_ = x_km_i_(2);
  double ct = cos(theta_km_);
  double st = sin(theta_km_);
  C_km_i_ << ct, st, -st, ct;

  /* Odometry 
  double t;
  input_k.get(u_k_km_, t);
  dp_k_km_ = u_k_km_.head(2);
  dtheta_k_km_ = u_k_km_(2);
  ct = cos(dtheta_k_km_);
  st = sin(dtheta_k_km_);
  C_k_km_ << ct, st, -st, ct;

  /* Step forward 
  p_k_i_ = p_km_i_ + C_km_i_.transpose() * dp_k_km_;
  C_k_i_ = C_k_km_ * C_km_i_;

  /* Write state at k 
  theta_k_ = acos( C_k_i_(0, 0) );
  x_k_i_.head(2) = p_k_i_;
  x_k_i_(2) = theta_k_;
  s_k.set(x_k_i_);*/
  
}

