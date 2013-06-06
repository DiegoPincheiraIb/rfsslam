// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include <Eigen/Core>
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
  Landmark(){};

  /** 
   * Constructor - defined only for our convenience and non-essential
   */
  Landmark(StateType x, UncertaintyType Sx);

  /** Default destructor */
  ~Landmark(){};

};
template<class StateType, class UncertaintyType> 
Landmark<StateType, UncertaintyType>::Landmark(StateType x,UncertaintyType Sx){
set(x,Sx);
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
