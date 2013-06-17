// Pose classes for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include <Eigen/Core>
#include <Eigen/LU>
#include "State.hpp"

/**
 * \class Pose
 * \brief Generic Pose class for storing vehicle pose information
 * \author Keith Leung
 */
template <class VecType>
class Pose : public State< VecType >
{
public:
  
  /** Default constructor */
  Pose() : t_(-1) {}

  /** Constructor */
  Pose(VecType &x, double t = -1) : t_(t) {
    State<VecType>::set(x);
  }

  /** Destructor */
  ~Pose(){}

  /** Get pose with time */
  void get( VecType &x, double &t){
    State<VecType>::get(x);
    t = t_;
  }

  /** Get pose without time */
  void get( VecType &x){
    State<VecType>::get(x);
  }

  /** Set pose with time */
  void set( VecType &x, double &t){
    State<VecType>::set(x);
    t_ = t;
  }

  /** Set pose without time */
  void set( VecType &x ){
    State<VecType>::set(x);
  }

protected:

  double t_; /**< time */
};

/**
 * \class PoseWithUncertainty
 * \brief Generic Pose class for storing vehicle pose information 
 *  with uncertainty
 * \author Keith Leung
 */
  template <class VecType, class MatType>
  class PoseWithUncertainty : public StateWithUncertainty< VecType, MatType >
{
public:
  
  /** Default constructor */
  PoseWithUncertainty() : t_(-1) {}

  /** Constructor */
  PoseWithUncertainty(VecType &x, double t = -1) : t_(t) {
    set(x, (MatType() << MatType::Zero()).finished() );
  }

  /** Constructor */
  PoseWithUncertainty(VecType &x, MatType &Sx, double t = -1) : t_(t) {
    set(x, Sx);
  }

  /** Destructor */
  ~PoseWithUncertainty(){}

  /** Get pose with time */
  void get( VecType &x, double &t){
    State<VecType>::get(x);
    t = t_;
  }

  /** Set pose with time */
  void set( VecType &x, double &t){
    State<VecType>::set(x);
    t_ = t;
  }

  /** Get pose */
  void get( VecType &x){
    State<VecType>::get(x);
  }

  /** Set pose */
  void set( VecType &x){
    State<VecType>::set(x);
  }

  /** Get pose and uncertainty with time */
  void get( VecType &x, MatType &Sx, double &t){
    StateWithUncertainty<VecType, MatType>::get(x, Sx);
    t = t_;
  }

  /** Get pose and uncertainty without time */
  void get( VecType &x, MatType &Sx){
    StateWithUncertainty<VecType, MatType>::get(x, Sx);
  }

  /** Set pose and uncertainty with time */
  void set( VecType &x, MatType &Sx, double &t){
    StateWithUncertainty<VecType, MatType>::set(x, Sx);
    t_ = t;
  }

  /** Set pose and uncertainty without time */
  void set( VecType &x, MatType &Sx ){
    StateWithUncertainty<VecType, MatType>::set(x, Sx);
  }

protected:

  double t_; /**< time */
  
};


/********** Define a 1d vehicle pose state **********/

/** Definition for 1d pose */
typedef PoseWithUncertainty< Eigen::Matrix<double, 1, 1>,
			     Eigen::Matrix<double, 1, 1> > Pose1d;

/********** Example implementation of a 2d vehicle pose state **********/

/**
 * \class Pose2d
 * \brief 2d vehicle Pose for (x,y) coordinate and rotation
 * \author Keith Leung
 * \note This class is derived from pose only so that we can add a 
 *  custom custructor for our convenience.
 */
class Pose2d : public PoseWithUncertainty< Eigen::Vector3d, Eigen::Matrix3d >
{

public:

  /** 
   * Default constructor, implementation of which can be empty
   */
  Pose2d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x pose sate
   * \param Sx pose state uncertainty (covariance matrix)
   * \param t time
   */
  Pose2d(Vec &x, double t = -1);

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x pose sate
   * \param Sx pose state uncertainty (covariance matrix)
   * \param t time
   */
  Pose2d(Vec &x, Mat &Sx, double t = -1);

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
