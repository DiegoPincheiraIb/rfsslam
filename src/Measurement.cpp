#include "Measurement.hpp"

#define PI 3.14159265359

/********** Implementation of example 1d measurement **********/

Measurement1d::Measurement1d(){}

Measurement1d::Measurement1d(double z, double Sz,double t)
{
  set(z,Sz,t);
}
Measurement1d::~Measurement1d(){}

double Measurement1d::mahalanobisDist2(double &z){

  double e = z - this->x_;
  if( this->pSxInv_ == NULL ){
    this->SxInv_ = 1 / this->Sx_;
    pSxInv_ = &(this->SxInv_);
  }

  return e * this->SxInv_ *e;
}

double Measurement1d::evaluateLikelihood(double &z){

  double md2 = mahalanobisDist2(z);
  return exp( -0.5 * md2 ) / ( sqrt( 2 * PI *Sx_) );
}

/********** Implementation of example 2d measurement **********/

Measurement2d::Measurement2d(){}

Measurement2d::~Measurement2d(){}

double Measurement2d::mahalanobisDist2(Eigen::Vector2d &z){

  Eigen::Vector2d e = z - this->x_;
  if(this->pSxInv_ == NULL){
    this->SxInv_ = this->Sx_.inverse();
  }
  return e.transpose() * this->SxInv_ * e;
}

double Measurement2d::evaluateLikelihood(Eigen::Vector2d &z){
  double mahalanobis_squared = mahalanobisDist2(z);
  return exp(-0.5 * mahalanobis_squared ) / (2 * PI * sqrt( Sx_.determinant() ) );

}





/********** Implementation of example 2d odometry measurement **********/

Odometry2d::Odometry2d(){}

Odometry2d::Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km, 
		       double vardx_k_km, double vardy_k_km, 
		       double vartheta_k_km, double t){
  Eigen::Vector3d u;
  Eigen::Vector3d covu_diag; 
  Eigen::Matrix3d covu;
  u << dx_k_km, dy_k_km, dtheta_k_km;
  covu_diag << vardx_k_km, vardy_k_km, vartheta_k_km;
  covu = covu_diag.asDiagonal();

  set( u, covu, t);
}

Odometry2d::~Odometry2d(){}









