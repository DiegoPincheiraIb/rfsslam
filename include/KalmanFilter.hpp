/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Class for using the Kalman filter equations
// Felipe Inostroza 2013

#ifndef KALMANFILTER_HPP
#define KALMANFILTER_HPP

#include "MeasurementModel.hpp"
#include "ProcessModel.hpp"

namespace rfs
{

/** 
 * \class KalmanFilter
 * \brief An Extended Kalman Filter (EKF) for performing estimate predictions and 
 * updates on a Landmark, given a MeasurementModel, ProcessModel, and a sensor Pose
 * \note In the future, we will augment this class with two more correct and update
 * functions so this can be used to update sensor pose, but there is no immediate
 * need for this.
 * \tparam ProcessModelType A process model derived from 
 * ProcessModel or StaticProcessModel
 * \tparam MeasurementModelType A measurement model derived from MeasurementModel
 * \author Felipe Inostroza, Keith Leung
 */
template <class ProcessModelType, class MeasurementModelType>
class KalmanFilter
{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef typename ProcessModelType::TInput TInput;
  typedef typename MeasurementModelType::TPose TPose;
  typedef typename MeasurementModelType::TLandmark TLandmark;
  typedef typename MeasurementModelType::TMeasurement TMeasurement;
  typedef MeasurementModelType TMeasurementModel;
  typedef ProcessModelType TProcessModel;

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
   * Kalman filter prediction step
   * This funtion uses the ProcessModel to propagate the Landmark estimate forward in time
   * \param[in] landmark_current The landmark before the prediction
   * \param[out] landmark_updated The landmark after the prediction step (this can be the same as landmark_current to update in place)
   * \param[in] dT time step size, if required by motion model
   */
  void predict(TLandmark &landmark_current, TLandmark &landmark_updated,
	       double dT = 1);

  /**
   * Kalman filter update step.
   * \note This function may be reimplemented in a derived class
   * \param[in] pose the sensor pose
   * \param[in] measurement the measurement to use for updating the landmark state
   * \param[in] landmark_current the current landmark state
   * \param[out] landmark_updated the updated landmark state
   * \param[out] zLikelihood if supplied, this stores the measurement likelihood 
   * (required by the RBPHDFilter)
   * \param[out] mahalanobisDist2 if supplied, stores the sqaured mahalanobis distance 
   * used to calculate zlikelihood (required by the RBPHDFilter)
   * \return true if the correction was performed. This could be false if calculateInnovation returns false
   */
  virtual bool correct(const TPose &pose, const TMeasurement &measurement, 
		       TLandmark &landmark_current, TLandmark &landmark_updated,
		       double* zLikelihood = NULL, 
		       double* mahalanobisDist2 = NULL);
  
  /**
   * Calculate the innovation. This may be reimplemented in a derived class
   * for special cases such as dealing with rotations.
   * \param[in] z_exp expected measurement 
   * \param[in] z_act actual measurement
   * \param[out] z_innov innovation
   * \return true to allow the correction step to proceed. This is useful for ignoring an update when
   * the innovation is too large due to outliers, which may cause problems for the filter.
   */
  virtual bool calculateInnovation(::Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_exp, 
				   ::Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_act,
				   ::Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_innov);

protected:

  MeasurementModelType *pMeasurementModel_;
  ProcessModelType *pProcessModel_;

  TMeasurement measurement_exp; /**< expected measurement */
  ::Eigen::Matrix <double, TLandmark::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime>  K; /**< Kalman gain */
  ::Eigen::Matrix <double, TMeasurement::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime > H; /**< measurement model Jacobian */
  ::Eigen::Matrix <double, TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> S; /**< innovation covariance */
  ::Eigen::Matrix <double, TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> S_inv; /**< inverse innovation covariance */
  ::Eigen::Matrix <double, TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime>  I; /**< identity matrix */
  ::Eigen::Matrix <double, TMeasurement::Vec::RowsAtCompileTime, 1> z_act; /**< measurement - actual */
  ::Eigen::Matrix <double, TMeasurement::Vec::RowsAtCompileTime, 1> z_exp; /**< measurement - expected */
  ::Eigen::Matrix <double, TMeasurement::Vec::RowsAtCompileTime, 1> z_innov; /**< measurement - innovation */
  ::Eigen::Matrix <double, TLandmark::Vec::RowsAtCompileTime, 1> m; /**< mean */
  ::Eigen::Matrix <double, TLandmark::Vec::RowsAtCompileTime, 1> m_updated; /**< mean - updated */
  ::Eigen::Matrix <double, TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime> P; /**< covariance */
  ::Eigen::Matrix <double, TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime> P_updated; /**< covariance - updated */
  RandomVec< ::Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, 1>, ::Eigen::Matrix < double , TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> > innov; /**< RandomVec form of innovation */

};

//********** Implementation of the standard Kalman Filter  **********/
// Felipe Inostroza 2013, Keith Leung 2014

template <class ProcessModelType, class MeasurementModelType> 
KalmanFilter<ProcessModelType, MeasurementModelType>::
KalmanFilter() 
{
  pMeasurementModel_ = NULL;
  pProcessModel_ = NULL;
}

template <class ProcessModelType, class MeasurementModelType> 
KalmanFilter<ProcessModelType, MeasurementModelType>::
KalmanFilter(ProcessModelType *pProcessModel, 
	     MeasurementModelType *pMeasurementModel)
{
  pMeasurementModel_= pMeasurementModel;
  pProcessModel_= pProcessModel;
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
predict(TLandmark &landmark_current, 
	TLandmark &landmark_updated,
	double dT){

  pProcessModel_->staticStep(landmark_updated, landmark_current, dT);

}

template <class ProcessModelType, class MeasurementModelType> 
bool KalmanFilter<ProcessModelType, MeasurementModelType>::
correct(const TPose &pose, const TMeasurement &measurement, 
	TLandmark &landmark_current, TLandmark &landmark_updated,
	double* zLikelihood, double* mahalanobisDist2 ){

  
  if(!pMeasurementModel_->measure( pose , landmark_current , measurement_exp , &H))
    return false; // invalid expected measurement produced
  
  TimeStamp t; 
  measurement.get(z_act, t);
  measurement_exp.get(z_exp, S);

  if(!calculateInnovation(z_exp, z_act, z_innov))
    return false; // innovation cannot be used, abort update
  
  landmark_current.get(m, P);

  I.setIdentity();
  measurement_exp.getCovInv(S_inv);
  K = P * H.transpose() * S_inv;
  P_updated = ( I - K * H ) * P;
  P_updated = ( P_updated + P_updated.transpose() ) / 2; // make sure covariance is symetric
  
  m_updated = m + K * z_innov;
  landmark_updated.set(m_updated, P_updated);
  
  if(zLikelihood != NULL){
    innov.set(z_act, S);
    *zLikelihood = innov.evalGaussianLikelihood( z_exp, mahalanobisDist2 );   
    if(*zLikelihood != *zLikelihood) // When likelihood is so small that it becomes NAN
      *zLikelihood = 0;  
  }

  return true;

}


template <class ProcessModelType, class MeasurementModelType> 
bool KalmanFilter<ProcessModelType, MeasurementModelType>::
calculateInnovation(::Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_exp, 
		    ::Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_act,
		    ::Eigen::Matrix< double, TMeasurement::Vec::RowsAtCompileTime, 1> &z_innov){
  z_innov = z_act - z_exp;
  return true;
}

}
#endif
