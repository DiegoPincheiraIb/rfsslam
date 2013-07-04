// Pose classes for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include "State.hpp"

/********** Define a 1d vehicle pose state **********/

/** Definition for 1d pose */
typedef RandomVec< Eigen::Matrix<double, 1, 1>,
		   Eigen::Matrix<double, 1, 1> > Pose1d;

/********** Example implementation of a 2d vehicle pose state **********/

/**
 * \class Pose2d
 * \brief 2d vehicle Pose for (x,y) coordinate and rotation
 * \author Keith Leung
 * \note This class is derived from pose only so that we can add a 
 *  custom custructor for our convenience.
 */
class Pose2d : public RandomVec< Eigen::Vector3d, Eigen::Matrix3d >
{

public:

  /** 
   * Default constructor, implementation of which can be empty
   */
  Pose2d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x pose vector
   * \param Sx pose uncertainty covariance
   * \param t time
   */
  Pose2d(Vec &x, Mat &Sx, double t = -1);

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x pose vector
   * \param Sx pose uncertainty covariance 
   * \param t time
   */
  Pose2d(Vec &x, double t = -1);

  /**
   * Constructor - defined only for our convenience and non-essential
   * \param x x-position
   * \param y y-position
   * \param theta orientation
   * \param va_x x-position variance
   * \param va_y y-position variance
   * \param va_theta theta orientation variance
   * \param t time
   */
  Pose2d( double x, double y, double theta, 
	  double var_x = 0, double var_y = 0, double var_theta = 0,
	  double t = -1 ); 

  /** Default destructor */
  ~Pose2d();

};


#endif
