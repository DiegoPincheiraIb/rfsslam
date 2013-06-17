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
 * \brief An abstract class for measurements.
 * \tparam VecType Eigen vector representing a measurement
 * \tparam MatType Eigen matrix representing 
 *  measurement uncertainties
 * \author Keith Leung, Felipe Inostroza
 */
template<class VecType, class MatType>
class Measurement 
  : public StateWithUncertainty<VecType, MatType>
{
public:

  typedef VecType MeasureType;
  typedef MatType MeasureUncertaintyType; 

  /** 
   * Default constructor   
   * \param nDim number of dimensions in measurement vector
   */
  Measurement(){ 
    t_ = -1; 
  }

  /** 
   * Constructor
   * \param nDim number of dimensions in measurement vector
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time information. 
   */
  Measurement(VecType &z, 
	      MatType &Sz, double t=-1){
    set(z, Sz, t);
  }

  /** 
   * Constructor
   * \param nDim number of dimensions in measurement vector
   * \param z - measurement
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time information. 
   */
  Measurement(VecType &z, 
	      double t=-1){
    set(z, (MatType() << MatType::Zero()).finished() , t);
  }

  /** Default destructor */
  ~Measurement(){};

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time information.
   */
  void set(VecType &z, 
	   MatType &Sz, 
	   double t = -1)
  {
    StateWithUncertainty< VecType, 
			  MatType >::set(z, Sz);
    t_ = t;
  }

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param t - time at wich the measurement was taken, negative  times indicate absence of time imformation.
   */
  void set(VecType &z, double t = -1)
  {
    State< VecType >::set(z);
    t_ = t;
  }

  /** 
   * Function for getting measurement values
   * \param z - measurement [overwritten]
   * \param Sz - measurement uncertainty [overwritten]
   * \param t - time at wich the measurement was taken, 
   *            negative times indicate absence of time imformation. [overwritten]
   */
  void get(VecType &z, MatType &Sz, double &t)
  {
    StateWithUncertainty< VecType, MatType >::get(z, Sz);
    t = t_;
  }

  /** 
   * Function for getting the value of the measurement
   * \param z measurement [overwritten]
   * \param t timestamp of the measurement, negative times indicate absence of time imformation.[overwritten]
   */
  void get(VecType &z, double &t)
  {
    State< VecType >::get(z);
    t = t_;
  }

protected:

  double t_; /**< Timestamp for the measurement, negative values indicate absence of time information.*/

};


/********** Define a 1d measurement **********/

/** Definition for 1d measurement */
typedef Measurement < Eigen::Matrix<double, 1, 1>,
		      Eigen::Matrix<double, 1, 1> > Measurement1d;

/********** Define a 2d measurement **********/

/** Definition for 2d measurement */
typedef Measurement <Eigen::Vector2d, Eigen::Matrix2d> Measurement2d;


/********** Examples implementation of 2d odometry measurement **********/

/**
 * \class Odometry2d
 * \brief A class for 2d odometry measurements for a 2d motion model
 * \author Keith Leung
 * \note This class is derived from measurement so that we can make a customized constructor.
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
