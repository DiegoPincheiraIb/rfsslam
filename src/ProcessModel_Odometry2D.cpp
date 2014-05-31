/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "ProcessModel_Odometry2D.hpp"

using namespace rfs;

MotionModel_Odometry2d::MotionModel_Odometry2d(){}

MotionModel_Odometry2d::MotionModel_Odometry2d( Pose2d::Mat &Q ) : ProcessModel(Q) {}

MotionModel_Odometry2d::~MotionModel_Odometry2d(){}

void MotionModel_Odometry2d::step(  Pose2d &s_k, 
				   Pose2d &s_km, 
				   Odometry2d &input_k, 
				   TimeStamp const &dT ){
 

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
  TimeStamp t_km;
  s_km.get(x_km_i_, t_km);
  p_km_i_ = x_km_i_.head(2);
  theta_km_ = x_km_i_(2);
  double ct = cos(theta_km_);
  double st = sin(theta_km_);
  C_km_i_ << ct, st, -st, ct;

  /* Odometry */
  input_k.get(u_k_km_);
  dp_k_km_ = u_k_km_.head(2);
  dtheta_k_km_ = u_k_km_(2);
  ct = cos(dtheta_k_km_);
  st = sin(dtheta_k_km_);
  C_k_km_ << ct, st, -st, ct;

  /* Step forward */
  p_k_i_ = p_km_i_ + C_km_i_.transpose() * dp_k_km_;
  C_k_i_ = C_k_km_ * C_km_i_;

  /* Write state at k */
  TimeStamp t_k = t_km + dT;
  theta_k_ = atan2( C_k_i_(0, 1), C_k_i_(0, 0) );
  x_k_i_.head(2) = p_k_i_;
  x_k_i_(2) = theta_k_;
  s_k.set(x_k_i_, t_k);
}

