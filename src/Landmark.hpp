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

  /** Default destructor */
  ~Landmark(){};

};




class Landmark1d : public Landmark<double, double>
{
public: 
  Landmark1d();
  ~Landmark1d();
};

/********** Example implementation of a 2d Landmark **********/

/**
* \class Landmark2d
* \brief 2d Landmark with no signature
* \author Felipe Inostroza
*/
class Landmark2d : public Landmark<Eigen::Vector3d, Eigen::Matrix3d>
{
public:
  /** Default constructor */
  Landmark2d();


  /** 
   * Constructor - defined only for our convenience and non-essential
   */

  Landmark2d(Eigen::Vector3d,Eigen::Matrix3d);

  /** Default destructor */
  ~Landmark2d();

  
};

#endif
