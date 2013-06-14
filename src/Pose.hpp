// Pose classes for defining vehicle state
// Keith Leung 2013

#ifndef POSE_HPP
#define POSE_HPP

#include <Eigen/Core>
#include "State.hpp"

/**
 * \class Pose
 * \brief Generic Pose class for storing vehicle pose information
 * \author Keith Leung
 */
template <class StateType>
class Pose : public State< StateType >
{
public:
  
  /** Default constructor */
  Pose() : t_(-1) {}

  /** Constructor */
  Pose(StateType x, double t = -1) : t_(t) {
    State<StateType>::set(x);
  }

  /** Destructor */
  ~Pose(){}

  /** Get pose with time */
  void get( StateType &x, double &t){
    State<StateType>::get(x);
    t = t_;
  }

  /** Get pose without time */
  void get( StateType &x){
    State<StateType>::get(x);
  }

  /** Set pose with time */
  void set( StateType x, double &t){
    State<StateType>::set(x);
    t_ = t;
  }

  /** Set pose without time */
  void set( StateType x ){
    State<StateType>::set(x);
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
  template <class StateType, class UncertaintyType>
  class PoseWithUncertainty : public StateWithUncertainty< StateType, UncertaintyType >
{
public:
  
  /** Default constructor */
  PoseWithUncertainty() : t_(-1) {}

  /** Constructor */
  PoseWithUncertainty(StateType x, UncertaintyType Sx, double t = -1) : t_(t) {
    StateWithUncertainty<StateType, UncertaintyType>::set(x, Sx);
  }

  /** Destructor */
  ~PoseWithUncertainty(){}

  /** Get pose with time */
  void get( StateType &x, double &t){
    State<StateType>::get(x);
    t = t_;
  }

  /** Get pose without time */
  void get( StateType &x){
    State<StateType>::get(x);
  }

  /** Get pose and uncertainty with time */
  void get( StateType &x, UncertaintyType &Sx, double &t){
    StateWithUncertainty<StateType, UncertaintyType>::get(x, Sx);
    t = t_;
  }

  /** Get pose and uncertainty without time */
  void get( StateType &x, UncertaintyType &Sx){
    StateWithUncertainty<StateType, UncertaintyType>::get(x, Sx);
  }

  /** Set pose and uncertainty with time */
  void set( StateType x, UncertaintyType Sx, double &t){
    StateWithUncertainty<StateType, UncertaintyType>::set(x, Sx);
    t_ = t;
  }

  /** Set pose and uncertainty without time */
  void set( StateType x, UncertaintyType Sx ){
    StateWithUncertainty<StateType, UncertaintyType>::set(x, Sx);
  }

protected:

  double t_; /**< time */
  
};


/********** Define a 1d vehicle pose state **********/

/** Definition for 1d pose */
typedef Pose< Eigen::Matrix<double, 1, 1> > Pose1d;


/********** Example implementation of a 2d vehicle pose state **********/

/**
 * \class Pose2d
 * \brief 2d vehicle Pose for (x,y) coordinate and rotation
 * \author Keith Leung
 * \note This class is derived from pose only so that we can add a 
 *  custom custructor for our convenience.
 */
class Pose2d : public Pose<Eigen::Vector3d>
{

public:

  typedef Eigen::Vector3d StateType;

  /** 
   * Default constructor, implementation of which can be empty
   */
  Pose2d();

  /** 
   * Constructor - defined only for our convenience and non-essential
   * \param x pose
   * \param t time
   */
  Pose2d(StateType x, double t = -1);

  /**
   * Constructor - defined only for our convenience and non-essential
   * \param x x-position
   * \param y y-position
   * \param theta rotation
   * \param t time
   */
  Pose2d( double x, double y, double theta, double t = -1 ); 

  /** Default destructor */
  ~Pose2d();

};


#endif
