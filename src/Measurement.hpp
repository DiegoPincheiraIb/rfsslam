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
 * \todo make measurement dimension a variable required in constructor so that 
 * a generic likelihood calculation function can be used
 *
 * \todo look into moving evaluateLikelihood into base class StateWithUncertainty
 *
 * \todo if we assume users will always use Eigen, mahalanobisDist2 can also move into
 * base class StateWithUncertainty
 *
 * \todo if we assume users will always use Eigen, store the inverse uncertainty 
 * in the constructor and when set() is called
 */
template<class MeasurementValType, class MeasurementUncertaintyType>
class Measurement 
  : public StateWithUncertainty<MeasurementValType, MeasurementUncertaintyType>
{
public:

  typedef MeasurementValType MeasureType;
  typedef MeasurementUncertaintyType MeasureUncertaintyType; 

  /** 
   * Default constructor   
   * \param nDim number of dimensions in measurement vector
   */
  Measurement(unsigned int nDim = 0) 
    : StateWithUncertainty<MeasurementValType, MeasurementUncertaintyType>(nDim) { 
    t_ = -1; 
  }

  /** 
   * Constructor
   * \param nDim number of dimensions in measurement vector
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time information. 
   */
  Measurement(unsigned int nDim, MeasurementValType z, 
	      MeasurementUncertaintyType Sz, double t=-1)
    : StateWithUncertainty<MeasurementValType, MeasurementUncertaintyType>( nDim )
  {
    set(z, Sz, t);
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
    StateWithUncertainty< MeasurementValType, MeasurementUncertaintyType >::set(z, Sz);
    t_ = t;
  }

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time imformation.
   */
  void set(MeasurementValType z, double t = -1)
  {
    State< MeasurementValType >::set(z);
    t_ = t;
  }

  /** 
   * Function for getting measurement values
   * \param z - measurement [overwritten]
   * \param Sz - measurement uncertainty [overwritten]
   * \param t - time at wich the measurement was taken, 
   *            negative times indicate absence of time imformation. [overwritten]
   */
  void get(MeasurementValType &z, MeasurementUncertaintyType &Sz, double &t)
  {
    StateWithUncertainty< MeasurementValType, MeasurementUncertaintyType >::get(z, Sz);
    t = t_;
  }

  /** 
   * Function for getting the value of the measurement
   * \param z measurement [overwritten]
   * \param t timestamp of the measurement, negative times indicate absence of time imformation.[overwritten]
   */
  void get(MeasurementValType &z, double &t)
  {
    State< MeasurementValType >::get(z);
    t = t_;
  }

protected:

  double t_; /**< Timestamp for the measurement, negative values indicate absence of time information.*/

};


/********** Example implementation of a 1d measurement **********/

/**
 * \class Measurement1d
 * \brief 1d Measurement
 * \author Felipe Inostroza
 */
class Measurement1d 
  : public Measurement < Eigen::Matrix<double, 1, 1>,
			 Eigen::Matrix<double, 1, 1> >
{
public:

  /** 
   * Default constructor, implementation of which can be empty 
   */
  Measurement1d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Measurement1d( Eigen::Matrix<double, 1, 1> z, 
		 Eigen::Matrix<double, 1, 1> Sz, 
		 double t=-1);

  /** Default destructor */
  ~Measurement1d();

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
