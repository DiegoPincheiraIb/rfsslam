#include "Pose.hpp"

Pose2d::Pose2d(){
  x_ << 0, 0, 0;
}

Pose2d::Pose2d(vec x){
  setPose(x);
}

void Pose2d::setPose(vec x){
  x_ = x;
}

void Pose2d::getPose(vec& x){
  x = x_;
}
