// Measurement Classes
// Felipe Inostroza, Keith Leung  2013

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include "RandomVec.hpp"

/** Definition for a  NULL measurement */
typedef RandomVec < Eigen::Matrix<double, 1, 1>,
		      Eigen::Matrix<double, 1, 1> > NullInput;

/** Definition for 1d measurement */
typedef RandomVec < Eigen::Matrix<double, 1, 1>,
		      Eigen::Matrix<double, 1, 1> > Measurement1d;

/** Definition for 2d measurement */
typedef RandomVec <Eigen::Vector2d, Eigen::Matrix2d> Measurement2d;


/********** Examples implementation of 2d odometry measurement **********/

/**
 * \class Odometry2d
 * \brief A class for 2d odometry measurements for a 2d motion model
 * \author Keith Leung
 * \note This class is derived from so that we can make a customized constructor.
 */
class Odometry2d : public RandomVec< Eigen::Vector3d, Eigen::Matrix3d >
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  /** Default constructor */
  Odometry2d();
 
  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x measurement vector
   * \param Sx measurement uncertainty covariance
   * \param t time
   */
  Odometry2d(Vec &x, Mat &Sx, double t = -1);

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x measurement vector
   * \param Sx measurement uncertainty covariance 
   * \param t time
   */
  Odometry2d(Vec &x, double t = -1);


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
