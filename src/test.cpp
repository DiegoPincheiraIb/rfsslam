#include "GaussianMixture.hpp"
#include "ParticleFilter.hpp"
#include "MeasurementModel.hpp"
#include <iostream>


int main(int argc, char* argv[]){
  
  RangeBearingModel measurementModel;
  Eigen::Matrix2d cov;
  cov << 1,0,0,3;
  measurementModel.setCov(cov);
  Measurement<int, int> m;
  
  Particle<Pose2d> p;

  int nParticles = 10;
  Pose2d x_0(2, 1, 0);
  OdometryMotionModel2d motionModel;

  ParticleFilter<OdometryMotionModel2d, RangeBearingModel> pf(nParticles, x_0, &motionModel, &measurementModel); 

  GaussianMixture<Landmark2d> map;


  return 0;
}
