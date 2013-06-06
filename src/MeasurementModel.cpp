#include <math.h>
#include "MeasurementModel.hpp"




/********** Implementation of example 2d measurement Model (Range and Bearing) **********/

RangeBearingModel::RangeBearingModel(){

}


RangeBearingModel::~RangeBearingModel(){

}



void RangeBearingModel::getCov(Eigen::Matrix2d &covZ){

  
  covZ=covZ_;
}


void RangeBearingModel::setCov(Eigen::Matrix2d covZ){
  covZ_=covZ;
}

void RangeBearingModel::setCov(double Sr,double Sb){
  covZ_<< Sr,0,
          0 ,Sb;
}

void RangeBearingModel::predict(Pose2d  pose, Landmark2d landmark, Measurement2d &prediction){

  Eigen::Vector3d robotPose;
  Eigen::Vector2d mean,landmarkState;
  Eigen::Matrix2d H,landmarkUncertainty;
  double range, bearing;


  pose.get(robotPose);
  landmark.get(landmarkState,landmarkUncertainty);

  range= sqrt(pow(landmarkState(0)-robotPose(0),2)+pow(landmarkState(1)-robotPose(1),2));
  
  bearing=atan2(landmarkState(1)-robotPose(1),landmarkState(0)-robotPose(0))-robotPose(2);

  while(bearing>PI) bearing-=2*PI;
  while(bearing<-PI) bearing+=2*PI;

  mean << range , bearing ;

  H << -(landmarkState(1)-robotPose(1))/(pow(mean(0),2)) , (landmarkState(0)-robotPose(0))/pow(mean(0),2),
       (landmarkState(0)-robotPose(0))/mean(0)           , (landmarkState(1)-robotPose(1))/mean(0) ;
  
  
  prediction.set(mean,H*landmarkUncertainty*H.transpose()+covZ_);

}



void RangeBearingModel::inversePredict(Pose2d pose, Landmark2d &landmark,Measurement2d measurement){

  Eigen::Vector3d poseState;
  Eigen::Vector2d measurementState,mean;
  Eigen::Matrix2d measurementUncertainty,covariance,Hinv;
  double t;

  pose.get(poseState);
  measurement.get(measurementState,measurementUncertainty,t);
  mean << poseState(0)+measurementState(0)*cos(poseState(2)+measurementState(1)),
          poseState(1)+measurementState(0)*sin(poseState(2)+measurementState(1));

  Hinv << cos(poseState(2)+measurementState(1)) , -measurementState(0)*sin(poseState(2)+measurementState(1)) ,
          sin(poseState(2)+measurementState(1)) , measurementState(0)*cos(poseState(2)+measurementState(1));


  landmark.set(mean,Hinv*measurementUncertainty*Hinv.transpose());

}


double RangeBearingModel::evaluateLikelyhood(Measurement2d prediction,Measurement2d measurement){

  Eigen::Vector2d predictionMean,measurementMean,diff;
  Eigen::Matrix2d predictionCov,measurementCov;
  double retval,t1,t2,det,mahalanobis;
  prediction.get(predictionMean,predictionCov,t1);
  measurement.get(measurementMean, measurementCov,t2);
  diff=measurementMean-predictionMean;
  det=predictionCov.determinant();
  mahalanobis=diff.transpose()*predictionCov.inverse()*diff;
  retval=exp(-0.5*mahalanobis)/(2*PI*sqrt(det));
  return retval;



}



