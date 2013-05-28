// Pose class for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include <Eigen/Core>

/**
 * \class Pose
 * \brief An abstract base class for defining the vehicle pose state
 * \author Keith Leung
 */
class Pose
{

public:

  /** Default constructor */
  Pose();

  /** Default destructor */
  ~Pose();

  /** Abstract function for setting the pose state */
  virtual void setPose();
  
  /** Abstract function for getting the pose state */
  virtual void getPose();

private:

};

/**
 * \class Pose2d
 * \brief 2d vehicle pose for (x,y) coordinate and rotation
 * \author Keith Leung
 */
class Pose2d : public Pose
{

public:

  typedef Eigen::Vector3d vec;

  /** Default constructor */
  Pose2d();

  /** Constructor */
  Pose2d(vec x);

  /** Default destructor */
  ~Pose2d();

  /** 
   * Function for setting pose state 
   * \param pose the pose to set
   */
  void setPose(vec pose);

  /** 
   * Function for getting pose state 
   * \param pose variable passed in by reference to be overwritten
   */
  void getPose(vec& pose);

private:

  vec x_;

};

#endif
