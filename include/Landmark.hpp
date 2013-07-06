// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include "RandomVec.hpp"

/** 
 * \class Landmark
 * \brief An abstract class for defining landmark state
 * \author Keith Leung
 */
template<class VecType, class MatType>
class Landmark : public RandomVec<VecType, MatType>
{
public:

  /** Default constructor */
  Landmark(){ 
    nReferences_ = 0; 
  }

  /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Landmark(VecType x, MatType Sx){
    this->set(x, Sx);
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

typedef Landmark< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >
Landmark1d;

typedef Landmark<Eigen::Vector2d, Eigen::Matrix2d> Landmark2d;

#endif
