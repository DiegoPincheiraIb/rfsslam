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

  typedef PoseType TPose;
  typedef LandmarkType TLandmark;
  typedef MeasurementType TMeasurement;
  
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
   * \param jacobian Jacobian of the measurement model at the point where the prediction is made [overwritten]
   * \todo there should be a overloaded version of this without the Jacobian
   */
  virtual void predict( PoseType &pose, LandmarkType &landmark, 
			MeasurementType &prediction , Eigen::Matrix<double , MeasurementType::Vec::RowsAtCompileTime ,
		  LandmarkType::Vec::RowsAtCompileTime > &jacobian ) = 0;

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

  /**
   * Abstract function of determining a landmark's probability of detection
   * \note If this is not reimplemented in a derived class, it will always
   * return a probability of detection of 1
   * \brief In addition to specifying a single probability of detection as the
   * return value, we can also specify whether a landmark is close to the sensing limit,
   * where there can be mismatch between the probability of detection returned by this
   * function, and the actual probability of detection in real life. This flag is read 
   * by the RB-PHD-Filter and is important in its implementation to prevent landmarks
   * near the sensing limit from dissapearing due to the probability of detection mismatch.
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isClostToSensingLimit true if landmark is close to the sensing limit,
   * where there can be a mismatch between the returned value and the actual probability
   * of detection in reality.
   * \return probability of detection
   */
  virtual double probabilityOfDetection( PoseType &pose,
					 LandmarkType &landmark,
					 bool &isCloseToSensingLimit){
    isCloseToSensingLimit = false;
    return 1;
  }


  /**
   * Abstract function for Determining the clutter intensity in the measurement space
   * \note This should be reimplemented in a derived class
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  virtual double clutterIntensity( MeasurementType &z,
				   int nZ ){
    return 0;
  }

  /**
   * Abstract function for Determining the clutter intensity integral
   * \note This should be reimplemented in a derived class
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity integral
   */
  double clutterIntensityIntegral( int nZ ){
    return 0;
  }  

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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for the 2d RangeBearingModel */
  struct Config{
    double probabilityOfDetection_;
    double probabilityOfFalseAlarm_; /**<  interpret as p( NULL | measurement exists) */
    double rangeLim_; /**< sensing range limit, beyond which Pd = 0 */
    double rangeLimBuffer_; /**< A buffer for Pd ambiguity */
  }config;

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
   * \param jacobian  jacobian of the measurement model evaluated at the prediction
   */
  void predict(Pose2d &pose, Landmark2d &landmark, Measurement2d &prediction, Eigen::Matrix2d &jacobian);

  /** 
   * Inverse measurement model 
   * \param pose robot pose 
   * \param landmark predicted landmark location [overwritten]
   * \param measurement measurement
   */
  void inversePredict(Pose2d &pose, Measurement2d &measurement, Landmark2d &landmark);

  /** 
   * Determine the probability of detection
   * \note This is where we can indirectly specify sensing limits and other sensor characteristics
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] Pd_upper
   * \param[out] Pd_lower
   * \return probability of detection
   */
  double probabilityOfDetection( Pose2d &pose,
				 Landmark2d &landmark,
				 bool &isCloseToSensingLimit);

  /**
   * Determine the clutter intensity in measurement space
   * \note uniform clutter intensity is assumed
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  double clutterIntensity( Measurement2d &z,
			   int nZ );

  /**
   * Determine the clutter intensity integral in measurement space
   * \brief This is calculated based on the probablity of false alarm,
   * defined as p( NULL | measurement exists)
   * \param[in] nZ the cardinality of Z
   * \return clutter intensity
   */
  double clutterIntensityIntegral( int nZ );

  
private:

  Eigen::Matrix2d covZ_;
  double sensingArea_;

};

#endif
