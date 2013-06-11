// Class for containing measurements with uncertainty
// Felipe Inostroza 2013

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP


#include <Eigen/LU>
#include <Eigen/Core>
#include <math.h>
#include "State.hpp"

/** 
 * \class Measurement
 * \brief An abstract class for measurements
 * \tparam MeasurementValType Object type representing measurements
 * \tparam MeasurementUncertaintyType Object type representing measurement uncertainties
 * \author Keith Leung, Felipe Inostroza
 *
 */
template<class MeasurementValType, class MeasurementUncertaintyType>
class Measurement : public StateWithUncertainty<MeasurementValType,MeasurementUncertaintyType>
{
public:

  /** 
   * Default constructor   
   */
  Measurement(){ t_ = -1; }

  /** 
   * Constructor
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time information. 
   */
  Measurement(MeasurementValType z, MeasurementUncertaintyType Sz, double t=-1)
  {
   
    t_=t;
    this->x_=z;
    this->Sx_=Sz;
  }

  /** Default destructor */
  ~Measurement(){};

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time information.
   */
  void set(MeasurementValType z, MeasurementUncertaintyType Sz, double t = -1)
  {
    this->x_ = z;
    this->Sx_ = Sz;
    t_ = t;
  }

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time imformation.
   */
  void set(MeasurementValType z, double t = -1)
  {
    this->x_ = z;
    t_ = t;
  }

  /** 
   * Function for getting measurement values
   * \param z - measurement [overwritten]
   * \param Sz - measurement uncertainty [overwritten]
   * \param t - time at wich the measurement was taken, 
   *            negative times indicate absence of time imformation. [overwritten]
   */
  void get(MeasurementValType &z, MeasurementUncertaintyType &Sz, double &t){
    z = this->x_;
    Sz = this->Sx_;
    t = t_;
  }

  /** 
   * Function for getting the value of the measurement
   * \param z measurement [overwritten]
   * \param t timestamp of the measurement, negative times indicate absence of time imformation.[overwritten]
   */
  void get(MeasurementValType &z, double &t){
    z = this->x_;
    t = t_;
  }

  /** 
   * Abstract function for returning the Mahalanobis distance from this measurement
   * \param z the measurement to which we measure the distance to
   * \return mahalanobis distance
   */
  virtual double mahalanobisDist( MeasurementValType &z){
    return -1;
  };


  /** 
   * Function for returning the Mahalanobis distance from this measurement
   * \param z the measurement to which we measure the distance to
   * \return mahalanobis distance
   */
  double mahalanobisDist(Measurement &z){
    return mahalanobisDist(z.x_);
  };

  /** 
   * Abstract function for returning the likelihood of a measurement
   * \param z the measurement whose likelihood will be evaluated
   * \return likelihood
   */
  virtual double evaluateLikelihood(MeasurementValType &z){
    return -1;
  };

  /** 
   * Function for returning the likelihood of a measurement
   * \param z the measurement whose likelihood will be evaluated
   * \return likelihood
   */
  double evaluateLikelihood(Measurement &z){
    return  evaluateLikelihood(z.x_);
  };


protected:

  double t_; /**< Timestamp for the measurement, negative values indicate absence of time information.*/

};


/********** Example implementation of a 1d measurement **********/

/**
 * \class Measurement1d
 * \brief 1d Measurement
 * \author Felipe Inostroza
 */
class Measurement1d : public Measurement <double,double>
{
public:

  /** 
   * Default constructor, implementation of which can be empty 
   */
  Measurement1d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Measurement1d(double x, double Sx, double t=-1);

  /** Default destructor */
  ~Measurement1d();

  /** 
   * Function for returning the Mahalanobis distance from this measurement
   * \param z the measurement to which we measure the distance to
   * \return mahalanobis distance
   */
  double mahalanobisDist( double &z);

  /** 
   * Function for returning the likelihood of a measurement
   * \param z the measurement whose likelihood will be evaluated
   * \return likelihood
   */
  double evaluateLikelihood(double &z);
};


/********** Example implementation of a 2d measurement **********/

/**
 * \class Measurement2d
 * \brief 2d Measurement
 * \author Felipe Inostroza
 */
class Measurement2d : public Measurement <Eigen::Vector2d, Eigen::Matrix2d>
{
public:
  /** 
   * Default constructor, implementation of which can be empty 
   */
  Measurement2d();

  /** Default destructor */
  ~Measurement2d();

  /** 
   * Function for returning the Mahalanobis distance from this measurement
   * \param z the measurement to which we measure the distance to
   * \return mahalanobis distance
   *
   * \todo include option for enforcing measurement range limit
   */
  double mahalanobisDist(Eigen::Vector2d  &z);

  /** 
   * Function for returning the likelihood of a measurement
   * \param z the measurement whose likelihood will be evaluated
   * \return likelihood
   */
  double evaluateLikelihood(Eigen::Vector2d &z);
};


/********** Examples implementation of 2d odometry measurement **********/

/**
 * \class Odometry2d
 * \brief A class for 2d odometry measurements for a 2d motion model
 * \author Keith Leung
 */
class Odometry2d : public Measurement< Eigen::Vector3d, Eigen::Matrix3d >
{
public:
  
  /** Default constructor */
  Odometry2d();
  
  /** 
   * Constructor, not necessary, but defined for convenience
   * \param dx_k_km x-displacement from frame k-1
   * \param dy_k_km y-displacement from frame k-1
   * \param dtheta_k_km rotational displacement from frame k-1
   * \param vardx_k_km variance in dx_k_km
   * \param vardy_k_km variance in dy_k_km
   * \param vardtheta_k_km variance in dtheta_k_km
   * \param t time of odometry reading
   */
  Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km,
	     double vardx_k_km, double vardy_k_km, double vartheta_k_km,
	     double t);

  /** Destructor */
  ~Odometry2d();
};

#endif
