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
 *     * Neither the name of the owner nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "GaussianMixture.hpp"
#include "ParticleFilter.hpp"
#include "MeasurementModel.hpp"
#include "RBPHDFilter.hpp"
#include <iostream>


int main(int argc, char* argv[]){
  
  RangeBearingModel measurementModel;
  Eigen::Matrix2d cov;
  cov << 1,0,0,3;
  measurementModel.setNoise(cov);

  
  Particle<Pose2d> p;

  int nParticles = 10;
  Pose2d x_0(2, 1, 0);
  Pose2d x_1;
  Odometry2d::Vec u;
  u << 0, 0, 0;
  Odometry2d odo;
  odo.set(u);

  Pose2d::Mat ProcessNoise;
   
  ProcessNoise << 3, 2, 1, 2, 4, -1, 1, -1, 5;

  OdometryMotionModel2d motionModel(ProcessNoise);
  motionModel.sample(x_1, x_0, odo);

  //ParticleFilter<OdometryMotionModel2d, RangeBearingModel> pf(nParticles, x_0, &motionModel, &measurementModel); 

  GaussianMixture<Landmark2d> map;

  // RBPHDFilter<OdometryMotionModel2d, RangeBearingModel> phdFilter(nParticles, x_0, &motionModel, &measurementModel);


  return 0;
}
