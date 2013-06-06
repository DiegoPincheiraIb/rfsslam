// Class to contain measurements
// Felipe Inostroza 2013

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include <Eigen/Core>

/** 
 * \class Measurement
 * \brief An abstract class for measurements
 * \tparam MeasurementValType Object type representing measurements
 * \tparam MeasurementUncertaintyType Object type representing measurement uncertainties
 * \author Keith Leung
 */

template<class MeasurementValType, class MeasurementUncertaintyType>
class Measurement
{
public:

  /** Default constructor */
  Measurement(){ t_ = 0; };


 /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Measurement(MeasurementValType z, MeasurementUncertaintyType Sz, double t=-1)
  {
    if(t >= 0)
      t_ = t;
    else
      t_=0;
    t_=t;
    z_=z;
    Sz_=Sz;
  };
  /** Default destructor */
  ~Measurement(){};

  /** 
   * Function for setting measurement values
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative or zero times indicate absence of time imformation.
   */
  void set(MeasurementValType z, MeasurementUncertaintyType Sz, double t = -1)
  {
    z_ = z;
    Sz_ = Sz;
    if(t >= 0)
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
    if(t >= 0)
      t_ = t;
  }

  /** 
   * Function for getting measurement values
   * \param z - measurement
   * \param Sz - measurement uncertainty
   * \param t - time at wich the measurement was taken, negative or zero times indicate absence of time imformation.
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

  Measurement1d(double x,double Sx,double t=-1);

  /** Default destructor */
  ~Measurement1d();

};


/********** Example implementation of a 2d measurement **********/

/**
 * \class Measurement2d
 * \brief 2d Measurement
 * \author Felipe Inostroza
 */
class Measurement2d : public Measurement <Eigen::Vector2d,Eigen::Matrix2d>
{
public:
  /** 
   * Default constructor, implementation of which can be empty 
   */
  Measurement2d();

  /** Default destructor */
  ~Measurement2d();

};



#endif
