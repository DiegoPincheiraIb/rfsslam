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

#include "MeasurementModel_6D.hpp"

namespace rfs
{

MeasurementModel_6D::MeasurementModel_6D(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}


MeasurementModel_6D::MeasurementModel_6D(::Eigen::Matrix3d &covZ){

  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_6D::MeasurementModel_6D(double Sx, double Sy, double Sz){

  Eigen::Matrix3d covZ;
  covZ <<  Sx, 0,  0,
		  0,   Sy, 0,
		  0,   0,  Sz;
  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLimMax_ = 5;
  config.rangeLimMin_ = 0.3;
  config.rangeLimBuffer_ = 0.25;
}

MeasurementModel_6D::~MeasurementModel_6D(){}

bool MeasurementModel_6D::measure(const Pose6d &pose,
				      const Landmark3d &landmark,
				      Measurement3d &measurement,
				      Eigen::Matrix3d *jacobian_wrt_lmk,
				      Eigen::Matrix<double, 3, 7> *jacobian_wrt_pose){


  Eigen::Vector3d mean, landmarkState;
  Eigen::Matrix3d H_lmk, landmarkUncertainty, cov;

  Eigen::Vector3d robotPosition;


  pose.getPos(robotPosition);
  landmark.get(landmarkState);
  Eigen::Quaterniond robotQ(pose.getRot());
  H_lmk = robotQ.toRotationMatrix();


  mean= H_lmk * (landmarkState-robotPosition);



  Eigen::Matrix<double, 3, 7> H_robot;
  Eigen::Matrix<double, 3, 4> H_robotrotation;
  Eigen::Matrix<double, 7, 7> robotUncertainty;
  double  range;

  pose.getCov( robotUncertainty);
  landmark.get(landmarkState,landmarkUncertainty);

  range = mean.norm();

  H_robotrotation <<   0           , 2*mean(2)  , -2*mean(1), 0,
		  	  	  	  -2*mean(2)  , 0          ,  2*mean(0), 0,
					  2*mean(1)   , -2*mean(0) ,  0        , 0; //skew symmetric matrix

  H_robot.block<3,3>(0,0) = -H_lmk;
  H_robot.block<3,4>(0,3) = H_robotrotation;


  cov = H_lmk * landmarkUncertainty * H_lmk.transpose() + H_robot * robotUncertainty * H_robot.transpose() + R_;
  measurement.set(mean, cov);

  if(jacobian_wrt_lmk != NULL)
    *jacobian_wrt_lmk = H_lmk;

  if(jacobian_wrt_pose != NULL)
    *jacobian_wrt_pose = H_robot;

  if(range > config.rangeLimMax_ || range < config.rangeLimMin_)
    return false;
  else
    return true;
}

void MeasurementModel_6D::inverseMeasure(const Pose6d &pose,
					 const Measurement3d &measurement,
					 Landmark3d &landmark){

  Eigen::Vector3d measurementState, mean;
  Eigen::Matrix3d measurementUncertainty, covariance, Hinv;

  Eigen::Vector3d robotPosition;

  pose.getPos(robotPosition);
  Eigen::Quaterniond robotQ(pose.getRot());
  Hinv = robotQ.conjugate().toRotationMatrix();


  measurement.get(measurementState);
  this->getNoise(measurementUncertainty);

  mean = Hinv*measurementState+robotPosition;


  covariance = Hinv * measurementUncertainty * Hinv.transpose();
  landmark.set( mean, covariance );

}

double MeasurementModel_6D::probabilityOfDetection( const Pose6d &pose,
						    const Landmark3d &landmark,
						    bool &isCloseToSensingLimit ){

  Pose6d::PosVec robotPose;
  Landmark3d::Vec landmarkState,diff;
  double range, Pd;

  isCloseToSensingLimit = false;

  pose.getPos(robotPose);
  landmark.get(landmarkState);
  diff=landmarkState-robotPose;

  range = diff.norm();

  if( range <= config.rangeLimMax_ && range >= config.rangeLimMin_){
    Pd = config.probabilityOfDetection_;
    if( range >= (config.rangeLimMax_ - config.rangeLimBuffer_ ) ||
	range <= (config.rangeLimMin_ + config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  }else{
    Pd = 0;
    if( range <= (config.rangeLimMax_ + config.rangeLimBuffer_ ) &&
	range >= (config.rangeLimMin_ - config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  }

  return Pd;
}

double MeasurementModel_6D::clutterIntensity( Measurement3d &z,
					    int nZ){
  return config.uniformClutterIntensity_;
}


double MeasurementModel_6D::clutterIntensityIntegral( int nZ ){
  double sensingVolume_ = 4.0/3.0 * PI * (pow(config.rangeLimMax_,3) - pow(config.rangeLimMin_,3));
  return ( config.uniformClutterIntensity_ * sensingVolume_ );
}

}
