// Measurement model
// Felipe Inostroza 2013

#ifndef MEASUREMENTMODEL_HPP
#define MEASUREMENTMODEL_HPP

#include <Eigen/Core>
#include "Measurement.hpp"
#include "Landmark.hpp"
#include "Pose.hpp"

/** 
 * \class MeasurementModel
 * \brief An abstract class for defining the measurement model
 * \author Felipe Inostroza, Keith Leung
 */
template<class PoseType, class LandmarkType, class MeasurementType>
class MeasurementModel
{
public:

  typedef PoseType tPose;
  typedef LandmarkType tLandmark;
  typedef MeasurementType tMeasurement;
  
  /** Default constructor */
  MeasurementModel(){};

  /** Default destructor */
  ~MeasurementModel(){};

  /** 
   * Abstract function for predicting measurements using pose and landmark estimates
   * \note This must be implemented in a derived class
   * \param pose robot pose 
   * \param landmark landmark wich will be used for predicting the measurement
   * \param prediction predicted measurement [overwritten]
   */
  virtual void predict( PoseType &pose, LandmarkType &landmark, 
			MeasurementType &prediction)=0;

 /** 
   * Abstract function for predicting landmark position from a robot pose and
   * its measurements
   * \note This must be implemented in a derived class
   * \param pose robot pose 
   * \param measurement measurement
   * \param landmark predicted landmark location [overwritten]
   */
  virtual void inversePredict( PoseType &pose,
			       MeasurementType &measurement, 
			       LandmarkType &landmark ) = 0;
};


////////// 2d Range-Bearing Measurement Model //////////

/** 
 * \class RangeBearingModel
 * \brief Range and bearing measurement model for 2d point landmarks.
 * The noise is considered to be Gaussian in the range and bearing space
 * \author Felipe Inostroza 
 */
                                                               
class RangeBearingModel : public MeasurementModel <Pose2d, Landmark2d, Measurement2d>{

public:

 /** Default constructor */
  RangeBearingModel();

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model
  * \param covZ measurement covariance
  */
  RangeBearingModel(Eigen::Matrix2d covZ);

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model, 
  * range and bearing are assumed to be uncorrelated
  * \param Sr Range variance
  * \param Sb Bearing variance
  */
  RangeBearingModel(double Sr, double Sb);

 /** Default destructor */
  ~RangeBearingModel();

 /**
  * Sets the uncertainty (covariance) of the measurement model
  * \param covZ measurement covariance
  */
  void setCov(Eigen::Matrix2d covZ);

 /**
  * Sets the uncertainty (covariance) of the measurement model, 
  * range and bearing are assumed to be uncorrelated
  * \param Sr Range variance
  * \param Sb Bearing variance
  */
  void setCov(double Sr, double Sb);

  /**
   * Function to get the covariance of the measurement model
   * \param covZ measurement covariance [overwritten]
   */
  void getCov(Eigen::Matrix2d &covZ);

  /** 
   * Predict a measurement from a pose and a landmark estimate
   * \param pose robot pose 
   * \param landmark landmark wich will be used for predicting the measurement
   * \param prediction predicted measurement [overwritten]
   */
  void predict(Pose2d &pose, Landmark2d &landmark, Measurement2d &prediction);

  /** 
   * Inverse measurement model 
   * \param pose robot pose 
   * \param landmark predicted landmark location [overwritten]
   * \param measurement measurement
   */
  void inversePredict(Pose2d &pose, Measurement2d &measurement, Landmark2d &landmark);
  
private:

  Eigen::Matrix2d covZ_;

};

#endif
