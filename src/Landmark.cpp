#include "Landmark.hpp"





/********** Implementation of example 1d vehicle pose state **********/
Landmark1d::Landmark1d(){
set(0,0);
}

Landmark1d::~Landmark1d(){}



/********** Implementation of example 2d vehicle pose state **********/


Landmark2d::Landmark2d(){
Eigen::Vector2d x;
x<<0,0;
Eigen::Matrix2d Sx;
Sx << 1,0,
     0,1;
set(x,Sx);

}

Landmark2d::~Landmark2d(){}







