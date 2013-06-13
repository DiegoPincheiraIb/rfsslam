// Pose classes for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include <Eigen/Core>
#include "State.hpp"

/********** Example implementation of a 2d vehicle pose state **********/

/**
 * \class Pose2d
 * \brief 2d vehicle Pose for (x,y) coordinate and rotation
 * \author Keith Leung
 */
class Pose2d : public State<Eigen::Vector3d>
{

public:

  typedef Eigen::Vector3d StateType;

  /** 
   * Default constructor, implementation of which can be empty
   */
  Pose2d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x - pose
   */
  Pose2d(StateType x);

  /**
   * Constructor - defined only for our convenience and non-essential
   * \param x x-position
   * \param y y-position
   * \param theta rotation
   */
  Pose2d( double x, double y, double theta ); 

  /** Default destructor */
  ~Pose2d();

};


/********** Example implementation of a 1d vehicle pose state **********/

/**
 * \class Pose1d
 * \brief 1d vehicle Pose
 * \author Keith Leung
 */
class Pose1d : public State<double>
{

public:

  /** 
   * Default constructor, implementation of which can be empty 
   */
  Pose1d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Pose1d(double x);

  /** Default destructor */
  ~Pose1d();

};

#endif
