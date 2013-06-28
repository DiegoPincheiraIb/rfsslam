// Class for using the Kalman filter equations
// Felipe Inostroza 2013

#ifndef KALMANFILTER_HPP
#define KALMANFILTER_HPP

#include "MeasurementModel.hpp"
#include "ProcessModel.hpp"
#include <Eigen/Core>
/** 
 * \class KalmanFilter
 * \brief A class for using the Kalman Filter Equations
 * \author Felipe Inostroza
 *
 * \todo Implement the Kalman Filter
 * 
 */
template <class PoseType, class LandmarkType, class MeasurementType> 
class KalmanFilter
{

public:

  /**
   * Default constructor
   */
  KalmanFilter();
  
  /**
   * Constructor
   * \param mM Pointer to the measurement model
   * \param pM Pointer to the process model
   */
  
  KalmanFilter(MeasurementModel< PoseType, LandmarkType, MeasurementType> *mM , ProcessModel< LandmarkType> *pM);
  
  /**
   * Function to set the the measurement model
   * \param mM Pointer to the measurement model to be used
   */
   
  void setMeasurementModel(MeasurementModel< PoseType, LandmarkType, MeasurementType> *mM);
  
  /**
   * Function to set the the process model
   * \param pM Pointer to the process model to be used
   */
   
  void setProcessModel(ProcessModel<LandmarkType> *pM);
  
  /**
   * Function to set the measurement and process models
   * 
   * \param mM Pointer to the measurement model
   * \param pM Pointer to the process model
   */
   
  void setModels(MeasurementModel< PoseType, LandmarkType, MeasurementType> *mM , ProcessModel< LandmarkType> *pM);
   

  /**
   * Function to apply the Kalman filter update step, this automatically checks whether it can 
   * reuse the Kalman gain or the predicted covariance.
   * \param new_landmark The updated landmark [overwritten]
   * \param prev_landmark The landmark that will be updated.
   * \param measurement The measurement used to update the state
   * \param model Measurement model that is used to calculate both the predicted measurement and its uncertainty
   *
   */
  
    

  
  virtual void correct(PoseType &pose , LandmarkType &new_landmark , LandmarkType &prev_landmark , MeasurementType &measurement);
  
  
   /**
   * This funtion uses the ProcessModel to propagate the feature through time, if no ProcessModel is set
   * nothing is done
   * \param new_landmark The landmark after the prediction step [overwritten]
   * \param prev_lanmark The landmark before the prediction  
   *
   */

  
  void predict(LandmarkType &new_landmark , LandmarkType &prev_landmark );
  
  /**
   * Function to calculate the innovation, when measuring rotations this needs to be overridden 
   * \param prediction Prediction of the measurement 
   * \param measurement measurement 
   *
   */
   
  virtual void calculateInnovation(Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime, 1>  &prediction , Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime, 1> &measurement);


protected:

MeasurementModel< PoseType, LandmarkType, MeasurementType> *measM_;
ProcessModel<LandmarkType> *procM_;
MeasurementType prediction_;

Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime , LandmarkType::Vec::RowsAtCompileTime>  K_ ;
Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime, LandmarkType::Vec::RowsAtCompileTime > H_ ;
Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime , LandmarkType::Vec::RowsAtCompileTime> S_ , Sinv_ , predCov_ , prevMeasCov_ ;
Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime , LandmarkType::Vec::RowsAtCompileTime>  newCov_ , prevLandmCov_ , I_;

Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime , 1> prevLandm_;
Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime, 1> pred_ , innovation_;
Eigen::Matrix < double , PoseType::Vec::RowsAtCompileTime , 1> prevPose_;


};

//********** Implementation of the standard Kalman Filter  **********/
// Felipe Inostroza 2013

template <class PoseType, class LandmarkType, class MeasurementType> 
KalmanFilter<PoseType , LandmarkType , MeasurementType>::
KalmanFilter(){
  measM_ = NULL;
  procM_ = NULL;
  I_.setIdentity();
}

template <class PoseType, class LandmarkType, class MeasurementType> 
KalmanFilter<PoseType , LandmarkType , MeasurementType>::
KalmanFilter(MeasurementModel< PoseType, LandmarkType, MeasurementType> *measM , ProcessModel<LandmarkType> *procM){
  measM_=measM;
  procM_=procM;
  I_.setIdentity();
}

template <class PoseType, class LandmarkType, class MeasurementType> 
void KalmanFilter<PoseType , LandmarkType , MeasurementType>::
setMeasurementModel(MeasurementModel< PoseType, LandmarkType, MeasurementType> *measM){

  measM_=measM;

}

template <class PoseType, class LandmarkType, class MeasurementType> 
void KalmanFilter<PoseType , LandmarkType , MeasurementType>::
setProcessModel(ProcessModel<LandmarkType> *procM){
  procM_=procM;
}

template <class PoseType, class LandmarkType, class MeasurementType> 
void KalmanFilter<PoseType , LandmarkType , MeasurementType>::
setModels(MeasurementModel< PoseType, LandmarkType, MeasurementType> *measM , ProcessModel<LandmarkType> *procM){

  measM_=measM;
  procM_=procM;
}

template <class PoseType, class LandmarkType, class MeasurementType> 
void KalmanFilter<PoseType , LandmarkType , MeasurementType>::
correct(PoseType &pose ,LandmarkType &new_landmark, LandmarkType &prev_landmark , MeasurementType &measurement){

  Eigen::Matrix < double , PoseType::Vec::RowsAtCompileTime , 1> newPose;
  Eigen::Matrix < double , MeasurementType::Vec::RowsAtCompileTime , 1> meas;
  Eigen::Matrix < double , MeasurementType::Vec::RowsAtCompileTime , MeasurementType::Vec::RowsAtCompileTime> measCov;
  Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime , 1> landm,newLandm;
  Eigen::Matrix < double , LandmarkType::Vec::RowsAtCompileTime , LandmarkType::Vec::RowsAtCompileTime> landmCov;
  
  double t;
  pose.get(newPose);
  prev_landmark.get(landm , landmCov);
  measurement.get(meas , measCov , t);
  
  if(newPose == prevPose_ && landm == prevLandm_ && landmCov == prevLandmCov_ ){
    
    if(measCov == prevMeasCov_){
    
      // Reuse K, newCov and measurement prediction
      
      
    }else{
    
      // Reuse the measurement prediction and its covariance, but not the innovation covariance (S) , 
      // gain (K) or new landmark Covariance (newCov) 
      
      S_ = predCov_+measCov;
      Sinv_ = S_.inverse();
      K_ = landmCov*H_.transpose()*Sinv_;
      newCov_ = ( I_ - K_*H_ ) * landmCov;
      
      prevMeasCov_=measCov;
    }
   
  }else{
    
    // Recalculate everything
    
    measM_->predict( pose , prev_landmark , prediction_ , H_ );
    prediction_.get(pred_ , predCov_ , t);
    S_ = predCov_+measCov;
    Sinv_ = S_.inverse();
    K_ = landmCov*H_.transpose()*Sinv_;
    newCov_ = ( I_ - K_*H_ ) * landmCov;
  

    prevPose_ = newPose;
    prevLandm_ = landm;
    prevLandmCov_ = landmCov;
  }
  
  calculateInnovation(pred_, meas);
  newLandm=landm+K_*innovation_;
  new_landmark.set(newLandm,newCov_);
  
  
}


template <class PoseType, class LandmarkType, class MeasurementType> 
void KalmanFilter<PoseType , LandmarkType , MeasurementType>::
calculateInnovation(Eigen::Matrix<double , LandmarkType::Vec::RowsAtCompileTime, 1> &prediction , Eigen::Matrix<double , LandmarkType::Vec::RowsAtCompileTime, 1> &measurement){


  innovation_=measurement-prediction;
  
  
  
}

template <class PoseType, class LandmarkType, class MeasurementType> 
void KalmanFilter<PoseType , LandmarkType , MeasurementType>::
predict(LandmarkType &new_landmark , LandmarkType &prev_landmark ){

  if(procM_!=NULL){
    procM_->step(new_landmark , prev_landmark);
  }else{
    new_landmark=prev_landmark;
  }

}




class RangeBearingKalmanFilter : public KalmanFilter <Pose2d , Landmark2d , Measurement2d>{

public:
  /**
   * Default constructor
   */
  RangeBearingKalmanFilter(){};
  
  /**
   * Constructor
   * \param mM Pointer to the measurement model
   * \param pM Pointer to the process model
   */
  
  RangeBearingKalmanFilter(MeasurementModel< Pose2d, Landmark2d, Measurement2d> *mM , ProcessModel< Landmark2d> *pM):
  KalmanFilter<Pose2d , Landmark2d , Measurement2d>(mM  , pM){};

  /**
   * Function to calculate the innovation 
   * \param prediction Prediction of the measurement 
   * \param measurement measurement 
   *
   */
   

  void calculateInnovation(Eigen::Vector2d &landmark , Eigen::Vector2d &measurement);

};

void RangeBearingKalmanFilter::calculateInnovation(Eigen::Vector2d &prediction , Eigen::Vector2d &measurement){



  innovation_=measurement-prediction;
  
  while(innovation_(1)>PI){
    innovation_(1)-=2*PI;
  }
  while(innovation_(1)<-PI){
    innovation_(1)+=2*PI;
  }  
  
    
}

typedef  KalmanFilter < Pose2d , Landmark2d , Measurement2d > KalmanFilter2d ;



#endif
