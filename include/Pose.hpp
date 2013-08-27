// Pose classes for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include "RandomVec.hpp"

/********** Define a 1d vehicle pose state **********/

/**
 * \class Pose1d
 * Vehicle position in 1d space
 * \brief Vehicle position
 * \author Keith Leung
 * \note The custom constructors are defined for convenience.
 */
class Pose1d : public RandomVec< Eigen::Matrix<double, 1, 1>,
				 Eigen::Matrix<double, 1, 1> >
{
public:
  
  /** For using Eigen fixed-size matrices with STL containers */
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** Default constructor */
  Pose1d();

  /** Constructor 
   *  \param[in] x position
   *  \param[in] Sx variance
   *  \param[in] t time
   */
  Pose1d(double x, double Sx, double t = -1);

  /** Constructor 
   *  \param[in] x position
   *  \param[in] Sx variance
   *  \param[in] t time
   */
  Pose1d(Eigen::Matrix<double, 1, 1> &x, Eigen::Matrix<double, 1, 1> &Sx, double t = -1);

  /** Constructor 
   *  \param[in] x position
   *  \param[in] t time
   */
  Pose1d(double x, double t = -1);

  /** Constructor 
   *  \param[in] x position
   *  \param[in] t time
   */
  Pose1d(Eigen::Matrix<double, 1, 1> &x, double t = -1);

  /** Destructor */
  ~Pose1d(){}

};

/********** 2d vehicle pose state **********/

/**
 * \class Pose2d
 * \f[ \mathbf{x} = \begin{bmatrix} x\\ y\\ \theta \end{bmatrix} \f]
 * \brief 2d vehicle Pose for (x,y) coordinate and rotation
 * \author Keith Leung
 * \note This class is derived from pose only so that we can add a 
 *  custom custructor for our convenience.
 */
class Pose2d : public RandomVec< Eigen::Vector3d, Eigen::Matrix3d >
{

public:

  /** For using Eigen fixed-size matrices with STL containers */
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** 
   * Default constructor
   */
  Pose2d();

  /** 
   * Constructor - defined only for our convenience
   * \param[in] x pose vector
   * \param[in] Sx pose uncertainty covariance
   * \param[in] t time
   */
  Pose2d(Vec &x, Mat &Sx, double t = -1);

  /** 
   * Constructor - defined only for our convenience
   * \param[in] x pose vector 
   * \param[in] t time
   */
  Pose2d(Vec &x, double t = -1);

  /**
   * Constructor - defined only for our convenience
   * \param[in] x x-position
   * \param[in] y y-position
   * \param[in] theta orientation
   * \param[in] va_x x-position variance
   * \param[in] va_y y-position variance
   * \param[in] va_theta theta orientation variance
   * \param[in] t time
   */
  Pose2d( double x, double y, double theta, 
	  double var_x = 0, double var_y = 0, double var_theta = 0,
	  double t = -1 ); 

  /** Default destructor */
  ~Pose2d();

};


#endif
