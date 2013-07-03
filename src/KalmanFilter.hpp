// Class for using the Kalman filter equations
// Felipe Inostroza 2013

#ifndef KALMANFILTER_HPP
#define KALMANFILTER_HPP

#include "MeasurementModel.hpp"
#include "ProcessModel.hpp"

/** 
 * \class KalmanFilter
 * \brief A Kalman Filter class for performing landmar estimate predictions and 
 * updates. It is intented for use in RBPHDFilter.
 * \tparam ProcessModelType process model derived from 
 * ProcessModel or StaticProcessModel
 * \tparam MeasurementModelType measurement model derived from MeasurementModel
 * \author Felipe Inostroza, Keith Leung
 */
template <class ProcessModelType, class MeasurementModelType>
class KalmanFilter
{

public:

  typedef typename ProcessModelType::TInput TInput;
  typedef typename MeasurementModelType::TPose TPose;
  typedef typename MeasurementModelType::TLandmark TLandmark;
  typedef typename MeasurementModelType::TMeasurement TMeasurement;

  /**
   * Default constructor
   */
  KalmanFilter();
  
  /**
   * Constructor
   * \param[in] pProcessModel Pointer to the process model
   * \param[in] pMeasurementModel Pointer to the measurement model
   */
  KalmanFilter(ProcessModelType *pProcessModel, 
	       MeasurementModelType *pMeasurementModel);
  
  /**
   * Set the the measurement model
   * \param[in] pMeasurementModel Pointer to the measurement model
   */
  void setMeasurementModel(MeasurementModelType *pMeasurementModel);
  
  /**
   * Set the the process model
   * \param[in] pProcessModel Pointer to the process model
   */
  void setProcessModel(ProcessModelType *pProcessModel);

  /**
   * Kalman filter update step. The function will  automatically check whether it can 
   * reuse the Kalman gain or the predicted covariance from a previous call.
   * \note This function may be reimplemented in a derived class
   * \param[in] pose The sensor pose
   * \param[in] measurement The measurement to use for updating the landmark state
   * \param[in]  landmark_current The current landmark state
   * \param[out] landmark_updated The updated landmark state
   * \param[out] zLikelihood if supplied, this stores the measurement likelihood 
   * (required by RBPHDFilter)
   * \warning this function is not thread-safe (if multiple threads using the same
   * instantiation of this class calls this function)
   */
  virtual void correct(TPose &pose, TMeasurement &measurement, 
		       TLandmark &landmark_current, TLandmark &landmark_updated,
		       double* zLikelihood = NULL);
  
   /**
   * This funtion uses the ProcessModel to propagate the feature through time
   * \param[in] lanmark_current The landmark before the prediction
   * \param[out] landmark_updated The landmark after the prediction step
   * \param[in] dT time step size, if required by motion model
   */
  void predict(TLandmark &landmark_current, TLandmark &landmark_updated,
	       double dT = 1);
  
  /**
   * Calculate the innovation. This may be reimplemented in a derived class
   * for special cases such as dealing with rotations. The result is stored in
   * innovation_
   * \param[in] z_exp Expected measurement 
   * \param[in] z_act Actual measurement
   */
  virtual void calculateInnovation(Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_exp, 
				   Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_act);


protected:

  MeasurementModelType *pMeasurementModel_;
  ProcessModelType *pProcessModel_;
  TMeasurement measurement_exp_;

  Eigen::Matrix < double , TLandmark::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime>  K_ ;
  Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime > H_ ;
  Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> S_ , S_inv_; 
  Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> z_exp_cov_ , R_prev_ ;
  Eigen::Matrix < double , TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime>  P_prev_ , I_;

  Eigen::Matrix < double , TLandmark::Vec::RowsAtCompileTime , 1> m_prev_;
  Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, 1> z_exp_ , innovation_;
  Eigen::Matrix < double , TPose::Vec::RowsAtCompileTime, 1> x_prev_;


};

//********** Implementation of the standard Kalman Filter  **********/
// Felipe Inostroza 2013

template <class ProcessModelType, class MeasurementModelType> 
KalmanFilter<ProcessModelType, MeasurementModelType>::
KalmanFilter() 
{
  pMeasurementModel_ = NULL;
  pProcessModel_ = NULL;
  I_.setIdentity();
}

template <class ProcessModelType, class MeasurementModelType> 
KalmanFilter<ProcessModelType, MeasurementModelType>::
KalmanFilter(ProcessModelType *pProcessModel, 
	     MeasurementModelType *pMeasurementModel)
{
  pMeasurementModel_= pMeasurementModel;
  pProcessModel_= pProcessModel;
  I_.setIdentity();
}

template <class ProcessModelType, class MeasurementModelType> 
void KalmanFilter<ProcessModelType, MeasurementModelType>::
setMeasurementModel(MeasurementModelType *pMeasurementModel){
  pMeasurementModel_ = pMeasurementModel;
}

template <class ProcessModelType, class MeasurementModelType> 
void KalmanFilter<ProcessModelType, MeasurementModelType>::
setProcessModel(ProcessModelType *pProcessModel){
  pProcessModel_ = pProcessModel;
}

template <class ProcessModelType, class MeasurementModelType> 
void KalmanFilter<ProcessModelType, MeasurementModelType>::
correct(TPose &pose, TMeasurement &measurement, 
	TLandmark &landmark_current, TLandmark &landmark_updated,
	double* zLikelihood){

  Eigen::Matrix < double , TPose::Vec::RowsAtCompileTime, 1> x;
  Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, 1> z_act;
  Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> R;
  Eigen::Matrix < double , TLandmark::Vec::RowsAtCompileTime, 1> m, m_updated;
  Eigen::Matrix < double , TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime> P, P_updated;;
  
  double t;
  pose.get( x );
  landmark_current.get(m , P);
  measurement.get(z_act , R , t);
  
  if(x == x_prev_ && m == m_prev_ && P == P_prev_ ){
    
    if(R == R_prev_){
    
      // Reuse K, newCov and measurement prediction
     
    }else{
    
      // Reuse the measurement prediction's covariance, 
      // but not the innovation covariance (S), 
      // gain (K), or landmark Covariance (P) 
      
      S_ = z_exp_cov_ + R;
      S_inv_ = S_.inverse();
      K_ = P * H_.transpose() * S_inv_;
      P_updated = ( I_ - K_*H_ ) * P;
      P_updated = ( P_updated + P_updated.transpose() ) / 2;
      
      R_prev_ = R;
    } 
   
  }else{
    
    // Recalculate everything
    
    pMeasurementModel_->measure( pose , landmark_current , measurement_exp_ , &H_ );
    measurement_exp_.get(z_exp_ , z_exp_cov_ , t);
    S_ = z_exp_cov_ + R;
    S_inv_ = S_.inverse();
    K_ = P * H_.transpose() * S_inv_;
    P_updated = ( I_ - K_*H_ ) * P;
    P_updated = ( P_updated + P_updated.transpose() ) / 2;
  
    x_prev_ = x;
    m_prev_ = m;
    P_prev_ = P;
  }
  
  calculateInnovation(z_exp_, z_act);
  m_updated= m + K_ * innovation_;
  landmark_updated.set(m_updated, P_updated);
  
  if(zLikelihood != NULL){
    StateWithUncertainty<Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, 1>, Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> > z_innov;
    z_innov.set( z_act, S_ );
    *zLikelihood = z_innov.evalGaussianLikelihood(z_exp_);
  }

}


template <class ProcessModelType, class MeasurementModelType> 
void KalmanFilter<ProcessModelType, MeasurementModelType>::
calculateInnovation(Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_exp, 
		    Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_act){
  innovation_ = z_act - z_exp;
}

template <class ProcessModelType, class MeasurementModelType> 
void KalmanFilter<ProcessModelType, MeasurementModelType>::
predict(TLandmark &landmark_current, 
	TLandmark &landmark_updated,
	double dT){

  pProcessModel_->staticStep(landmark_updated, landmark_current, dT);

}



/**
 * \class RangeBearingKalmanFilter
 * \brief A Kalman filter for updating a 2d landmark position from 
 * a single range-bearing measurements
 */
class RangeBearingKalmanFilter : 
  public KalmanFilter <StaticProcessModel<Landmark2d>, RangeBearingModel>{

  typedef RangeBearingModel::TMeasurement::Vec Vec;

public:
  /**
   * Default constructor
   */
  RangeBearingKalmanFilter(){};
  
  /**
   * Constructor
   * \param pMeasurementModel Pointer to the measurement model
   * \param pProcessModel Pointer to the process model
   */
  RangeBearingKalmanFilter(StaticProcessModel<Landmark2d> *pProcessModel,
			   RangeBearingModel *pMeasurementModel):
  KalmanFilter<StaticProcessModel<Landmark2d>, RangeBearingModel>
  (pProcessModel, pMeasurementModel){}

  /**
   * Function to calculate the innovation 
   * \param prediction Prediction of the measurement 
   * \param measurement measurement 
   */
  void calculateInnovation(Vec &z_exp, Vec &z_act){
    
    innovation_ = z_act - z_exp;
  
    while(innovation_(1)>PI){
      innovation_(1)-=2*PI;
    }
    while(innovation_(1)<-PI){
      innovation_(1)+=2*PI;
    }  
  }

};

#endif
