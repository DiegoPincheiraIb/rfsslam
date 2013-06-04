#include "Landmark.hpp"






/********** Implementation of example 1d vehicle pose state **********/
Landmark1d::Landmark1d(){}

Landmark1d::~Landmark1d(){}








/********** Implementation of example 2d vehicle pose state **********/


Landmark2d::Landmark2d(){
Eigen::Vector3d x;
x<<0,0,0;
Eigen::Matrix3d Sx;
Sx << 1,0,0,
     0,1,0,
     0,0,1;
set(x,Sx);

}

Landmark2d::~Landmark2d(){}

Landmark2d::Landmark2d(Eigen::Vector3d x,Eigen::Matrix3d Sx){
set(x,Sx);
}


