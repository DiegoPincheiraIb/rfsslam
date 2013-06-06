// Measurement model
// Felipe Inostroza 2013

#ifndef MEASUREMENTMODEL_HPP
#define MEASUREMENTMODEL_HPP



#include <Eigen/Core>
#include <Eigen/LU>
#include "Measurement.hpp"
#include "Landmark.hpp"
#include "Pose.hpp"
#define PI 3.14159265359

/** 
 * \class MeasurementModel
 * \brief An abstract class for defining the measurement model
 * \author Felipe Inostroza
 */
template<class PoseType,class LandmarkType, class MeasurementType>
class MeasurementModel
{
public:
 /** Default constructor */

  MeasurementModel(){};


 /** Default destructor */

  ~MeasurementModel(){};



 /** 
   * Abstract function for predicting measurements using both pose and landmark estimates
   * 
   * This must be implemented in a derived class
   * \param pose robot pose 
   * \param landmark landmark wich will be used for predicting the measurement
   * \param prediction predicted measurement [overwritten]
   */
  virtual void predict(PoseType pose,LandmarkType landmark,MeasurementType &prediction)=0;

 /** 
   * Abstract function for applying the inverse measurement model 
   * 
   * This must be implemented in a derived class
   * \param pose robot pose 
   * \param landmark predicted landmark location [overwritten]
   * \param measurement measurement
   */

  virtual void inversePredict(PoseType pose,LandmarkType &landmark,MeasurementType measurement)=0;

 /** 
   * Abstract function evaluating the likelyhood of a measruement 
   * 
   * This must be implemented in a derived class
   * \param prediction Measurement prediction, can be calculated using the predict function 
   * \param measurement measurement
   */


  virtual double evaluateLikelyhood(MeasurementType prediction,MeasurementType measurement)=0;

};
/** 
 * \class MeasurementModel
 * \brief Measurement model for a 2d landmark were the measurement the range and bearing to the object
 * (The noise is considered to be gaussian in the range and bearing space)
 * \author Felipe Inostroza
 */
                                                               
class RangeBearingModel : MeasurementModel <Pose2d,Landmark2d,Measurement2d>{

public:
 /** Default constructor */

  RangeBearingModel();


 /** Default destructor */

  ~RangeBearingModel();

 /**
  * Function to set the covariance of the measurement model
  * \param covZ measurement covariance
  *
  */
  void setCov(Eigen::Matrix2d covZ);
 /**
  * Function to set the covariance of the measurement model, range and bearing are assumed to be uncorrelated
  * \param Sr Range variance
  * \param Sb Bearing variance
  *
  */
  void setCov(double Sr,double Sb);


 /**
  * Function to get the covariance of the measurement model
  * \param covZ measurement covariance [overwritten]
  *
  */
  void getCov(Eigen::Matrix2d &covZ);



 /** 
   * Function for predicting measurements using both pose and landmark estimates
   * 
   * This must be implemented in a derived class
   * \param pose robot pose 
   * \param landmark landmark wich will be used for predicting the measurement
   * \param prediction predicted measurement [overwritten]
   */
  void predict(Pose2d pose,Landmark2d landmark,Measurement2d &prediction);

 /** 
   * Function for applying the inverse measurement model 
   * 
   * This must be implemented in a derived class
   * \param pose robot pose 
   * \param landmark predicted landmark location [overwritten]
   * \param measurement measurement
   */

  void inversePredict(Pose2d pose,Landmark2d &landmark,Measurement2d measurement);

 /** 
   * Function evaluating the likelyhood of a measruement 
   * 
   * This must be implemented in a derived class
   * \param prediction Measurement prediction, can be calculated using the predict function 
   * \param measurement measurement
   */


  double evaluateLikelyhood(Measurement2d prediction,Measurement2d measurement);



private:

Eigen::Matrix2d covZ_;


};

#endif
