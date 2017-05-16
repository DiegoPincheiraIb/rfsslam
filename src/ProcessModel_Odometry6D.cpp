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

#include "ProcessModel_Odometry6D.hpp"

using namespace rfs;

MotionModel_Odometry6d::MotionModel_Odometry6d(){}

MotionModel_Odometry6d::MotionModel_Odometry6d( Pose6d::Mat &Q ) : ProcessModel(Q) {}

MotionModel_Odometry6d::~MotionModel_Odometry6d(){}

void MotionModel_Odometry6d::step(  Pose6d &s_k,
				   Pose6d &s_km,
				   Odometry6d &input_k,
				   TimeStamp const &dT , Eigen::Matrix<double,7,7> *H){
 

  if (H != NULL){
    std::cerr << "MotionModel_Odometry6d: jacobian calculation not implemented!!!\n";
    exit(1);
  }
  Pose6d::PosVec p_k_i_;   /* \f[ \begin{bmatrix} x \\ y \\ z\end{bmatrix}_{k} \f]*/
  Pose6d::PosVec p_km_i_;  /* \f[ \begin{bmatrix} x \\ y \\ z\end{bmatrix}_{k-1} \f]*/
  Eigen::Quaterniond q_k_;          /* \f[ q_k \f\] */
  

  Eigen::Matrix3d Sd_k_km_; /* uncertainty of translation input */
  Odometry6d::PosVec dp_k_km_; /* translation input */

 
  /* State at k-1 */
  s_km.getPos(p_km_i_);
  Eigen::Quaterniond q_km_(s_km.getRot());          /* \f[ q_{k-1} \f\] */
  TimeStamp t_km=s_km.getTime();

  /* Odometry */
  input_k.getPos(dp_k_km_);
  Eigen::Quaterniond q_input_(input_k.getRot());  /* rotation Input */
  q_input_.normalize();

  /* Step forward */

  p_k_i_ = p_km_i_ + q_km_._transformVector(dp_k_km_);
  q_k_ =   q_km_ * q_input_;

  /* Write state at k */
  TimeStamp t_k = t_km + dT;




  s_k.setPos(p_k_i_);
  s_k.setRot(q_k_.coeffs(), t_k);
}

