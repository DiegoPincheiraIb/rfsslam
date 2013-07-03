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
