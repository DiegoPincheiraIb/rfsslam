// Measurement Classes
// Felipe Inostroza, Keith Leung  2013

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include "RandomVec.hpp"

/** Definition for a NULL measurement */
typedef RandomVec < Eigen::Matrix<double, 1, 1>,
		    Eigen::Matrix<double, 1, 1> > NullInput;

/** Definition for 1d measurement */
typedef RandomVec < Eigen::Matrix<double, 1, 1>,
		    Eigen::Matrix<double, 1, 1> > Measurement1d;

/** Definition for 2d measurement */
typedef RandomVec <Eigen::Vector2d, Eigen::Matrix2d> Measurement2d;


/** Definition for 1d odometry */
typedef RandomVec <  Eigen::Matrix<double, 1, 1>,
		     Eigen::Matrix<double, 1, 1> > Odometry1d;


/********** 2d odometry measurement **********/

/**
 * \class Odometry2d
 * \brief A class for 2d odometry measurements for a 2d motion model
 * \author Keith Leung
 * \note This class is derived from so that we can make a customized constructor for convenience.
 */
class Odometry2d : public RandomVec< Eigen::Vector3d, Eigen::Matrix3d >
{
public:

  /** For using Eigen fixed-size matrices with STL containers */
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  
  /** Default constructor */
  Odometry2d();
 
  /** 
   * Constructor - defined only for our convenience
   * \param[in] x measurement vector
   * \param[in] Sx measurement uncertainty covariance
   * \param[in] t time
   */
  Odometry2d(Vec &x, Mat &Sx, double t = -1);

  /** 
   * Constructor - defined only for our convenience
   * \param[in] x measurement vector
   * \param[in] t time
   */
  Odometry2d(Vec &x, double t = -1);


  /** 
   * Constructor - defined only for our convenience
   * \param[in] dx_k_km x-displacement from frame k-1
   * \param[in] dy_k_km y-displacement from frame k-1
   * \param[in] dtheta_k_km rotational displacement from frame k-1
   * \param[in] vardx_k_km variance in dx_k_km
   * \param[in] vardy_k_km variance in dy_k_km
   * \param[in] vardtheta_k_km variance in dtheta_k_km
   * \param[in] t time of odometry reading
   */
  Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km,
	     double vardx_k_km, double vardy_k_km, double vartheta_k_km,
	     double t);

  /** Destructor */
  ~Odometry2d();
};

#endif
