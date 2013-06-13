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

  /** Constructor
   *  \param nDim number of dimensions in landmark state
   */
  Landmark(unsigned int nDim = 0) 
    : StateWithUncertainty<StateType, UncertaintyType>(nDim){ 
    nReferences_ = 0; 
  }

  /** 
   * Constructor - defined only for our convenience and non-essential
   *  \param nDim number of dimensions in landmark state
   */
  Landmark(unsigned int nDim, StateType x, UncertaintyType Sx)
    : StateWithUncertainty<StateType, UncertaintyType>(nDim){
    set(x, Sx);
    nReferences_ = 0;
  }

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


protected:

  unsigned int nReferences_; /**< Number of references to this landmark */ 

};


/********** Example implementation of a 1d Landmark **********/

/**
* \class Landmark1d
* \brief 1d Landmark with no signature
* \author Felipe Inostroza
*/

class Landmark1d 
  : public Landmark< Eigen::Matrix<double, 1, 1>,
		     Eigen::Matrix<double, 1, 1> >
{
public: 
  /** Default constructor */
  Landmark1d();

  /** Default destructor */
  ~Landmark1d();

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

};








#endif
