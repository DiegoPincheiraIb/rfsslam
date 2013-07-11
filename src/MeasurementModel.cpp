#include <Eigen/LU>
#include <math.h>
#include "MeasurementModel.hpp"

/******** Implementation of example 2d measurement Model (Range and Bearing) **********/

RangeBearingModel::RangeBearingModel(){

  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLim_ = 5;
  config.rangeLimBuffer_ = 0.25;
}


RangeBearingModel::RangeBearingModel(Eigen::Matrix2d &covZ){

  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLim_ = 5;
  config.rangeLimBuffer_ = 0.25;
}

RangeBearingModel::RangeBearingModel(double Sr, double Sb){

  Eigen::Matrix2d covZ;
  covZ <<  Sr, 0, 0, Sb;
  setNoise(covZ);
  config.probabilityOfDetection_ = 0.95;
  config.uniformClutterIntensity_ = 0.1;
  config.rangeLim_ = 5;
  config.rangeLimBuffer_ = 0.25;
}

RangeBearingModel::~RangeBearingModel(){}

void RangeBearingModel::measure(Pose2d  &pose, Landmark2d &landmark, 
				Measurement2d &measurement, 
				Eigen::Matrix2d *jacobian){

  Eigen::Vector3d robotPose;
  Eigen::Vector2d mean, landmarkState;
  Eigen::Matrix2d H, landmarkUncertainty, cov;
  double range, bearing;

  pose.get(robotPose);
  landmark.get(landmarkState,landmarkUncertainty);

  range = sqrt(  pow(landmarkState(0) - robotPose(0), 2)
		+pow(landmarkState(1) - robotPose(1), 2) );
  bearing = atan2( landmarkState(1) - robotPose(1) , landmarkState(0) - robotPose(0) ) - robotPose(2);

  while(bearing>PI) bearing-=2*PI;
  while(bearing<-PI) bearing+=2*PI;

  mean << range, bearing ;

  H <<  (landmarkState(0)-robotPose(0))/mean(0)          , (landmarkState(1)-robotPose(1))/mean(0) ,
        -(landmarkState(1)-robotPose(1))/(pow(mean(0),2)), (landmarkState(0)-robotPose(0))/pow(mean(0),2) ;
  
  cov = H * landmarkUncertainty * H.transpose() ;
  measurement.set(mean, cov);
  
  if(jacobian != NULL)
    *jacobian = H;

}

void RangeBearingModel::inverseMeasure(Pose2d &pose, Measurement2d &measurement, 
				       Landmark2d &landmark){
  Eigen::Vector3d poseState;
  Eigen::Vector2d measurementState, mean;
  Eigen::Matrix2d measurementUncertainty, covariance, Hinv;
  double t;

  pose.get(poseState);
  measurement.get(measurementState, measurementUncertainty, t);
  mean << poseState(0) + measurementState(0) *cos( poseState(2) + measurementState(1) ),
          poseState(1) + measurementState(0) *sin( poseState(2) + measurementState(1) );

  Hinv << cos(poseState(2)+measurementState(1)) , -measurementState(0)*sin(poseState(2)+measurementState(1)) ,
          sin(poseState(2)+measurementState(1)) , measurementState(0)*cos(poseState(2)+measurementState(1));

  covariance = Hinv * measurementUncertainty *Hinv.transpose();
  landmark.set( mean, covariance );

}

double RangeBearingModel::probabilityOfDetection( Pose2d &pose,
						  Landmark2d &landmark,
						  bool &isCloseToSensingLimit ){

  Pose2d::Vec robotPose;
  Landmark2d::Vec landmarkState;
  double range, Pd;
  
  isCloseToSensingLimit = false;

  pose.get(robotPose);
  landmark.get(landmarkState);

  range = sqrt(  pow(landmarkState(0) - robotPose(0), 2)
		+pow(landmarkState(1) - robotPose(1), 2) );

  if( range <= config.rangeLim_){
    Pd = config.probabilityOfDetection_;
    if( range >= (config.rangeLim_ - config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  }else{
    Pd = 0;
    if( range <= (config.rangeLim_ + config.rangeLimBuffer_ ) )
      isCloseToSensingLimit = true;
  } 

  return Pd;
}

double RangeBearingModel::clutterIntensity( Measurement2d &z,
					    int nZ){
  return config.uniformClutterIntensity_;
}


double RangeBearingModel::clutterIntensityIntegral( int nZ ){
  double sensingArea_ = 2 * PI * config.rangeLim_;
  return ( config.uniformClutterIntensity_ * sensingArea_ );
}

