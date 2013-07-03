#include <math.h>
#include "ProcessModel.hpp"

/********** Implementation for an examples 2D motion model **********/

OdometryMotionModel2d::OdometryMotionModel2d(){}

OdometryMotionModel2d::OdometryMotionModel2d( Pose2d::Mat S ) : ProcessModel(S) {}

OdometryMotionModel2d::~OdometryMotionModel2d(){}

void OdometryMotionModel2d::step(  Pose2d &s_k, 
				   Pose2d &s_km, 
				   Odometry2d &input_k, 
				   const double dT ){
  
  /* State at k-1 */
  s_km.get(x_km_i_);
  p_km_i_ = x_km_i_.head(2);
  theta_km_ = x_km_i_(2);
  double ct = cos(theta_km_);
  double st = sin(theta_km_);
  C_km_i_ << ct, st, -st, ct;

  /* Odometry */
  double t;
  input_k.get(u_k_km_, t);
  dp_k_km_ = u_k_km_.head(2);
  dtheta_k_km_ = u_k_km_(2);
  ct = cos(dtheta_k_km_);
  st = sin(dtheta_k_km_);
  C_k_km_ << ct, st, -st, ct;

  /* Step forward */
  p_k_i_ = p_km_i_ + C_km_i_.transpose() * dp_k_km_;
  C_k_i_ = C_k_km_ * C_km_i_;

  /* Write state at k */
  theta_k_ = atan2( C_k_i_(0, 1), C_k_i_(0, 0) );
  x_k_i_.head(2) = p_k_i_;
  x_k_i_(2) = theta_k_;
  s_k.set(x_k_i_);
}

