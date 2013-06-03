#include "Pose.hpp"
#include "Landmark.hpp"
#include <iostream>


int main(int argc, char* argv[]){

  Pose2d pose;
  Pose2d::vec x, y;
  x << 1,2,3;
  pose.set(x);
  pose.get(y); 
  std::cout << y << std::endl;
  

  /*Landmark2d l;
  Landmark2d::vec lx;
  l.getState(lx);
  std::cout << "Landmark" << lx << std::endl;
  */
  return 0;
}
