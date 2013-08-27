#include <math.h>
#include "ProcessModel.hpp"

/********** Implementation for a 2D motion model **********/

OdometryMotionModel2d::OdometryMotionModel2d(){}

OdometryMotionModel2d::OdometryMotionModel2d( Pose2d::Mat &Q ) : ProcessModel(Q) {}

OdometryMotionModel2d::~OdometryMotionModel2d(){}

void OdometryMotionModel2d::step(  Pose2d &s_k, 
				   Pose2d &s_km, 
				   Odometry2d &input_k, 
				   const double dT ){
 

  Pose2d::Vec x_k_i_;       /* \f[ \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k} \f]*/
  Pose2d::Vec x_km_i_;      /* \f[ \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k-1} \f]*/
  Eigen::Vector2d p_k_i_;   /* \f[ \begin{bmatrix} x \\ y \end{bmatrix}_{k} \f]*/
  Eigen::Vector2d p_km_i_;  /* \f[ \begin{bmatrix} x \\ y \end{bmatrix}_{k-1} \f]*/
  double theta_k_;          /* \f[ \theta_k \f\] */
  double theta_km_;         /* \f[ \theta_k \f\] */
  Eigen::Matrix2d C_k_i_;   /* rotation matrix k */
  Eigen::Matrix2d C_km_i_;  /* rotation matrix k-1 */
  
  Eigen::Vector3d u_k_km_;  /* odometry input */
  Eigen::Matrix3d Sd_k_km_; /* uncertainty of translation input */
  Eigen::Vector2d dp_k_km_; /* translation input */
  double dtheta_k_km_;      /* rotation input */
  Eigen::Matrix2d C_k_km_;  /* rotation matrix from odometry input */
 
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


/************ Implementation of a 1D motion model ***********/

OdometryMotionModel1d::OdometryMotionModel1d( Pose1d::Mat &Q ) : ProcessModel(Q) {}

void OdometryMotionModel1d::step ( Pose1d &s_k, Pose1d &s_km, Odometry1d &input_k, 
				   double const dT){

  Pose1d::Vec x_km_; /**< position at k-1 */
  Pose1d::Vec x_k_;  /**< position at k */
  Odometry1d::Vec u_k_; /**< odometry from k-1 to k */

  /* k - 1 */
  s_km.get(x_km_);

  /* odometry */
  input_k.get(u_k_);

  /* step forward */
  x_k_ = x_km_ + u_k_;

  /* write state at k */
  s_k.set(x_k_);

}
