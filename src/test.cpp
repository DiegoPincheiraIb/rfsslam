#include "ParticleFilter.hpp"
#include <iostream>


int main(int argc, char* argv[]){

  Measurement<int, int> m;

  Particle<Pose2d> p;

  int nParticles = 10;
  Pose2d x_0(2, 1, 0);
  OdometryMotionModel2d motionModel;

  ParticleFilter<Pose2d, Odometry2d> pf(nParticles, x_0, &motionModel); 


  return 0;
}
