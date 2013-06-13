// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include <Eigen/Core>
#include <Eigen/LU>
#include "State.hpp"

/** 
 * \class Landmark
 * \brief An abstract class for defining landmark state
 * \author Keith Leung
 */
template<class StateType, class UncertaintyType>
class Landmark : public StateWithUncertainty<StateType, UncertaintyType>
{
public:

  /** Default constructor */
  Landmark(){ nReferences_ = 0; }

  /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Landmark(StateType x, UncertaintyType Sx);

  /** Default destructor */
  ~Landmark(){};

  /** 
   * Increase the count for objects referencing this object 
   * \return number of references 
   */
  unsigned int incNRef(){ nReferences_++; return getNRef(); }

  /**  
   * Decrease the count for objects referencing this object
   * \return number of references
   */
  unsigned int decNRef(){ nReferences_--; return getNRef(); }

  /**
   * Get the count of objects referencing this object
   * \return number of references
   */
  unsigned int getNRef(){ return nReferences_; }

  /** 
   * Abstract function for calculating the squared Mahalanobis distance from 
   * this object's state
   * \param x the state to which we measure the distance to
   * \return mahalanobis distance
   */
  virtual double mahalanobisDist2( StateType &x) = 0;

protected:

  unsigned int nReferences_; /**< Number of references to this landmark */ 

};

template<class StateType, class UncertaintyType> 
Landmark<StateType, UncertaintyType>::Landmark(StateType x,UncertaintyType Sx){
  set(x, Sx);
  nReferences_ = 0;
}



/********** Example implementation of a 1d Landmark **********/

/**
* \class Landmark1d
* \brief 1d Landmark with no signature
* \author Felipe Inostroza
*/

class Landmark1d : public Landmark<double, double>
{
public: 
 /** Default constructor */
  Landmark1d();

  /** Default destructor */
  ~Landmark1d();

   /** 
   *  Overwrite the Mahalanobis distance abstract (virtual) function
   *  \param x point to which we measure to
   *  \return the Mahalanobis distance
   */
   double mahalanobisDist2(double &x);
};

/********** Example implementation of a 2d Landmark **********/

/**
* \class Landmark2d
* \brief 2d Landmark with no signature
* \author Felipe Inostroza
*/
class Landmark2d : public Landmark<Eigen::Vector2d, Eigen::Matrix2d>
{
public:
  /** Default constructor */
  Landmark2d();

  /** Default destructor */
  ~Landmark2d();

  /** 
   *  Overwrite the Mahalanobis distance abstract (virtual) function
   *  \param x point to which we measure to
   *  \return the Mahalanobis distance
   */
  double mahalanobisDist2(Eigen::Vector2d &x);

};








#endif
