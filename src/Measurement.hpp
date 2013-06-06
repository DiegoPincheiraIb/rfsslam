// Class for containing measurements with uncertainty
// Felipe Inostroza 2013

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include <Eigen/Core>

/** 
 * \class Measurement
 * \brief An abstract class for measurements
 * \tparam MeasurementValType Object type representing measurements
 * \tparam MeasurementUncertaintyType Object type representing measurement uncertainties
 * \author Keith Leung, Felipe Inostroza
 *
 */
template<class MeasurementValType, class MeasurementUncertaintyType>
class Measurement
{
public:

  /** Default constructor 
   *  \note Should t_ default to 0 or -1? Theoretically we can have a measurement at 0? 
   */
  Measurement(){ t_ = 0; }

  /** 
   * Constructor
   * \todo add paramter list to comments
   * \todo if t < 0, it first gets set to 0, and then back to original value. Why?
   */
  Measurement(MeasurementValType z, MeasurementUncertaintyType Sz, double t=-1)
  {
    if(t >= 0)
      t_ = t;
    else
      t_ = 0;
    t_=t;
    z_=z;
    Sz_=Sz;
  }

  /** Default destructor */
  ~Measurement(){};

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative or zero times indicate absence of time information.
   */
  void set(MeasurementValType z, MeasurementUncertaintyType Sz, double t = -1)
  {
    z_ = z;
    Sz_ = Sz;
    t_ = t;
  }

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param t - time at wich the measurement was taken, negative or zero times indicate absence of time imformation.
   */
  void set(MeasurementValType u, double t = -1)
  {
    z_ = u;
    t_ = t;
  }

  /** 
   * Function for getting measurement values
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, 
   *            negative or zero times indicate absence of time imformation.
   */
  void get(MeasurementValType &z, MeasurementUncertaintyType &Sz, double &t){
    z = z_;
    Sz = Sz_;
    t = t_;
  }

  /** 
   * Function for getting process input values
   * \param u input [overwritten]
   * \param t time of input [overwritten]
   */
  void get(MeasurementValType &u, double &t){
    u = z_;
    t = t_;
  }


private:

  double t_; /**< time of the input */
  MeasurementValType z_;  /**< Input */
  MeasurementUncertaintyType Sz_; /**< Input covariance */

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
